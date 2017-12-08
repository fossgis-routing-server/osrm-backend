#include "engine/plugins/match.hpp"
#include "engine/plugins/plugin_base.hpp"

#include "engine/api/match_api.hpp"
#include "engine/api/match_parameters.hpp"
#include "engine/api/match_parameters_tidy.hpp"
#include "engine/map_matching/bayes_classifier.hpp"
#include "engine/map_matching/sub_matching.hpp"
#include "util/coordinate_calculation.hpp"
#include "util/integer_range.hpp"
#include "util/json_util.hpp"
#include "util/string_util.hpp"

#include <cstdlib>

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

namespace osrm
{
namespace engine
{
namespace plugins
{

// Filters PhantomNodes to obtain a set of viable candiates
void filterCandidates(const std::vector<util::Coordinate> &coordinates,
                      MatchPlugin::CandidateLists &candidates_lists)
{
    for (const auto current_coordinate : util::irange<std::size_t>(0, coordinates.size()))
    {
        bool allow_uturn = false;

        if (coordinates.size() - 1 > current_coordinate && 0 < current_coordinate)
        {
            double turn_angle =
                util::coordinate_calculation::computeAngle(coordinates[current_coordinate - 1],
                                                           coordinates[current_coordinate],
                                                           coordinates[current_coordinate + 1]);

            // sharp turns indicate a possible uturn
            if (turn_angle <= 90.0 || turn_angle >= 270.0)
            {
                allow_uturn = true;
            }
        }

        auto &candidates = candidates_lists[current_coordinate];
        if (candidates.empty())
        {
            continue;
        }

        // sort by forward id, then by reverse id and then by distance
        std::sort(candidates.begin(),
                  candidates.end(),
                  [](const PhantomNodeWithDistance &lhs, const PhantomNodeWithDistance &rhs) {
                      return lhs.phantom_node.forward_segment_id.id <
                                 rhs.phantom_node.forward_segment_id.id ||
                             (lhs.phantom_node.forward_segment_id.id ==
                                  rhs.phantom_node.forward_segment_id.id &&
                              (lhs.phantom_node.reverse_segment_id.id <
                                   rhs.phantom_node.reverse_segment_id.id ||
                               (lhs.phantom_node.reverse_segment_id.id ==
                                    rhs.phantom_node.reverse_segment_id.id &&
                                lhs.distance < rhs.distance)));
                  });

        auto new_end =
            std::unique(candidates.begin(),
                        candidates.end(),
                        [](const PhantomNodeWithDistance &lhs, const PhantomNodeWithDistance &rhs) {
                            return lhs.phantom_node.forward_segment_id.id ==
                                       rhs.phantom_node.forward_segment_id.id &&
                                   lhs.phantom_node.reverse_segment_id.id ==
                                       rhs.phantom_node.reverse_segment_id.id;
                        });
        candidates.resize(new_end - candidates.begin());

        if (!allow_uturn)
        {
            const auto compact_size = candidates.size();
            for (const auto i : util::irange<std::size_t>(0, compact_size))
            {
                // Split edge if it is bidirectional and append reverse direction to end of list
                if (candidates[i].phantom_node.forward_segment_id.enabled &&
                    candidates[i].phantom_node.reverse_segment_id.enabled)
                {
                    PhantomNode reverse_node(candidates[i].phantom_node);
                    reverse_node.forward_segment_id.enabled = false;
                    candidates.push_back(
                        PhantomNodeWithDistance{reverse_node, candidates[i].distance});

                    candidates[i].phantom_node.reverse_segment_id.enabled = false;
                }
            }
        }

        // sort by distance to make pruning effective
        std::sort(candidates.begin(),
                  candidates.end(),
                  [](const PhantomNodeWithDistance &lhs, const PhantomNodeWithDistance &rhs) {
                      return lhs.distance < rhs.distance;
                  });
    }
}

void ElongateInternalRouteResult(InternalRouteResult &base_result,
                                 InternalRouteResult &new_result,
                                 std::size_t ups_index)
{
    BOOST_ASSERT(base_result.is_valid());
    BOOST_ASSERT(new_result.is_valid());
    // In an InternalRouteResult object, each entry of the unpacked_path_segments vector
    // represents a leg of the route. When merging two InternalRouteResult objects, it's
    // ambiguous to which leg the new object's unpacked_path_segments should be added.
    // The `ups_index` value indcates to which unpacked_path_segment entry to elongate
    BOOST_ASSERT(ups_index <= base_result.unpacked_path_segments.size());
    // start a new leg in the base_result if the ups_index doesn't exist yet
    bool new_leg = ups_index == base_result.unpacked_path_segments.size();
    if (new_leg)
        base_result.unpacked_path_segments.push_back({});

    // collect all the iterators for writing
    auto new_segments_begin = new_result.unpacked_path_segments[ups_index].begin();
    auto new_src_reverse_begin = new_result.source_traversed_in_reverse.begin();
    auto new_trg_reverse_begin = new_result.target_traversed_in_reverse.begin();
    auto new_src_reverse_end = new_result.source_traversed_in_reverse.end();
    auto new_trg_reverse_end = new_result.target_traversed_in_reverse.end();
    auto new_segments_end = new_result.unpacked_path_segments[ups_index].end();
    // if the base result path segments end where the new ones' start, deduplicate by
    // moving the source write iterators (from the new_result) forward one
    if (base_result.unpacked_path_segments[ups_index].back().turn_via_node ==
        new_result.unpacked_path_segments[ups_index]
            .front()
            .turn_via_node) // FIXME there must be a better way to make this check
    {
        new_segments_begin++;
        new_src_reverse_begin++;
        new_trg_reverse_begin++;
    }
    std::move(new_segments_begin,
              new_segments_end,
              std::back_inserter(base_result.unpacked_path_segments[ups_index]));
    std::move(new_src_reverse_begin,
              new_src_reverse_end,
              std::back_inserter(base_result.source_traversed_in_reverse));
    std::move(new_trg_reverse_begin,
              new_trg_reverse_end,
              std::back_inserter(base_result.target_traversed_in_reverse));
    base_result.shortest_path_weight =
        base_result.shortest_path_weight + new_result.shortest_path_weight;
    if (new_leg)
    {
        // if this is a new leg, insert the phantom_node pair of the new_result into the
        // base_result segment_end_coordinates
        BOOST_ASSERT(base_result.segment_end_coordinates.back().target_phantom ==
                     new_result.segment_end_coordinates.front().source_phantom);
        BOOST_ASSERT(base_result.segment_end_coordinates.size() - 1 == ups_index);
        base_result.segment_end_coordinates.push_back(
            std::move(new_result.segment_end_coordinates.front()));
    }
    else
    {
        // otherwise, replace the last segment_end_coordinate with the last one in the new_result
        // segment_end_coordinates
        base_result.segment_end_coordinates.back().target_phantom =
            new_result.segment_end_coordinates.front().target_phantom;
    }
}

Status MatchPlugin::HandleRequest(const RoutingAlgorithmsInterface &algorithms,
                                  const api::MatchParameters &parameters,
                                  util::json::Object &json_result) const
{
    if (!algorithms.HasMapMatching())
    {
        return Error("NotImplemented",
                     "Map matching is not implemented for the chosen search algorithm.",
                     json_result);
    }

    if (!CheckAlgorithms(parameters, algorithms, json_result))
        return Status::Error;

    const auto &facade = algorithms.GetFacade();

    BOOST_ASSERT(parameters.IsValid());

    // enforce maximum number of locations for performance reasons
    if (max_locations_map_matching > 0 &&
        static_cast<int>(parameters.coordinates.size()) > max_locations_map_matching)
    {
        return Error("TooBig", "Too many trace coordinates", json_result);
    }

    if (!CheckAllCoordinates(parameters.coordinates))
    {
        return Error("InvalidValue", "Invalid coordinate value.", json_result);
    }

    // Check for same or increasing timestamps. Impl. note: Incontrast to `sort(first,
    // last, less_equal)` checking `greater` in reverse meets irreflexive requirements.
    const auto time_increases_monotonically = std::is_sorted(
        parameters.timestamps.rbegin(), parameters.timestamps.rend(), std::greater<>{});

    if (!time_increases_monotonically)
    {
        return Error(
            "InvalidValue", "Timestamps need to be monotonically increasing.", json_result);
    }

    SubMatchingList sub_matchings;
    api::tidy::Result tidied;
    if (parameters.tidy)
    {
        // Transparently tidy match parameters, do map matching on tidied parameters.
        // Then use the mapping to restore the original <-> tidied relationship.
        tidied = api::tidy::tidy(parameters);
    }
    else
    {
        tidied = api::tidy::keep_all(parameters);
    }

    // assuming radius is the standard deviation of a normal distribution
    // that models GPS noise (in this model), x3 should give us the correct
    // search radius with > 99% confidence
    std::vector<double> search_radiuses;
    if (tidied.parameters.radiuses.empty())
    {
        search_radiuses.resize(tidied.parameters.coordinates.size(),
                               routing_algorithms::DEFAULT_GPS_PRECISION * RADIUS_MULTIPLIER);
    }
    else
    {
        search_radiuses.resize(tidied.parameters.coordinates.size());
        std::transform(tidied.parameters.radiuses.begin(),
                       tidied.parameters.radiuses.end(),
                       search_radiuses.begin(),
                       [](const boost::optional<double> &maybe_radius) {
                           if (maybe_radius)
                           {
                               return *maybe_radius * RADIUS_MULTIPLIER;
                           }
                           else
                           {
                               return routing_algorithms::DEFAULT_GPS_PRECISION * RADIUS_MULTIPLIER;
                           }

                       });
    }

    auto candidates_lists = GetPhantomNodesInRange(facade, tidied.parameters, search_radiuses);

    filterCandidates(tidied.parameters.coordinates, candidates_lists);
    if (std::all_of(candidates_lists.begin(),
                    candidates_lists.end(),
                    [](const std::vector<PhantomNodeWithDistance> &candidates) {
                        return candidates.empty();
                    }))
    {
        return Error("NoSegment",
                     std::string("Could not find a matching segment for any coordinate."),
                     json_result);
    }

    // call the actual map matching
    sub_matchings =
        algorithms.MapMatching(candidates_lists,
                               tidied.parameters.coordinates,
                               tidied.parameters.timestamps,
                               tidied.parameters.radiuses,
                               parameters.gaps == api::MatchParameters::GapsType::Split);

    if (sub_matchings.size() == 0)
    {
        return Error("NoMatch", "Could not match the trace.", json_result);
    }

    // hardcode for now
    // TODO parse user supplied waypoints parameter
    std::vector<std::size_t> parameter_waypoints = {1, 5};
    for (const auto waypoint : parameter_waypoints)
    {
        bool found = false;
        std::for_each(sub_matchings.begin(), sub_matchings.end(), [&](const SubMatching &sm) {
            auto index = std::find(sm.indices.begin(), sm.indices.end(), waypoint);
            if (index != sm.indices.end())
                found = true;
        });
        if (!found)
            return Error(
                "NoMatch", "Requested waypoint parameter could not be matched.", json_result);
    }

    std::vector<InternalRouteResult> sub_routes(sub_matchings.size());
    // each sub_route will correspond to a MatchObject
    for (auto index : util::irange<std::size_t>(0UL, sub_matchings.size()))
    {
        BOOST_ASSERT(sub_matchings[index].nodes.size() > 1);

        // FIXME we only run this to obtain the geometry
        // The clean way would be to get this directly from the map matching plugin
        PhantomNodes current_phantom_node_pair;
        for (unsigned i = 0; i < sub_matchings[index].nodes.size() - 1; ++i)
        {
            // start a new leg if the sub_matching index being processed is a user requested
            // waypoint
            // TODO ... make this into a function
            bool new_leg = std::find(parameter_waypoints.begin(),
                                     parameter_waypoints.end(),
                                     sub_matchings[index].indices[i]) != parameter_waypoints.end();
            current_phantom_node_pair.source_phantom = sub_matchings[index].nodes[i];
            current_phantom_node_pair.target_phantom = sub_matchings[index].nodes[i + 1];
            BOOST_ASSERT(current_phantom_node_pair.source_phantom.IsValid());
            BOOST_ASSERT(current_phantom_node_pair.target_phantom.IsValid());
            // force uturns to be on
            // we split the phantom nodes anyway and only have bi-directional phantom nodes for
            // possible uturns
            std::vector<PhantomNodes> current_coords{current_phantom_node_pair};
            InternalRouteResult current_sub_result =
                algorithms.ShortestPathSearch(current_coords, {false});
            if (i == 0)
            {
                // start new route object, entry in sub_routes should be empty
                BOOST_ASSERT(!sub_routes[index].is_valid());
                sub_routes[index] = current_sub_result;
            }
            else if (new_leg)
            {
                BOOST_ASSERT(sub_routes[index].is_valid());
                // start new leg in InternalRouteResult
                auto number_of_legs = sub_routes[index].unpacked_path_segments.size();
                ElongateInternalRouteResult(
                    sub_routes[index], current_sub_result, number_of_legs + 1);
            }
            else
            {
                BOOST_ASSERT(sub_routes[index].is_valid());
                // merge new step into existing InternalRouteResult route leg
                auto number_of_legs = sub_routes[index].unpacked_path_segments.size();
                ElongateInternalRouteResult(sub_routes[index], current_sub_result, number_of_legs);
            }
        }
        BOOST_ASSERT(sub_routes[index].shortest_path_weight != INVALID_EDGE_WEIGHT);
    }

    api::MatchAPI match_api{facade, parameters, tidied};
    match_api.MakeResponse(sub_matchings, sub_routes, json_result);

    return Status::Ok;
}
}
}
}
