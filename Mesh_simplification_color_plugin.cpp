#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "ui_Mesh_simplification_color_dialog.h"

#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>

#include <boost/foreach.hpp>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

namespace SMS = CGAL::Surface_mesh_simplification;

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;

class Custom_stop_predicate {
	bool m_and;
	CGAL::Surface_mesh_simplification::Count_stop_predicate<FaceGraph> m_count_stop;
	CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double> m_length_stop;

public:
	Custom_stop_predicate(bool use_and, std::size_t nb_edges, double edge_length)
		: m_and(use_and),
		m_count_stop(nb_edges),
		m_length_stop(edge_length) {
		std::cout << "\nSimplifying until:" << std::endl
			<< " * Number of edges = " << nb_edges << std::endl;

		if (edge_length != std::numeric_limits<double>::max()) {
			std::cout << (use_and ? " AND " : " OR ") << std::endl
				<< " * Minimum edge length = " << edge_length << std::endl;
		}
	}

	template <typename Profile>
	bool operator() (const double& current_cost, const Profile& edge_profile,
		std::size_t initial_count, std::size_t current_count) const {
		if (m_and)
			return (m_count_stop(current_cost, edge_profile, initial_count, current_count)
				&& m_length_stop(current_cost, edge_profile, initial_count, current_count));
		else
			return (m_count_stop(current_cost, edge_profile, initial_count, current_count)
				|| m_length_stop(current_cost, edge_profile, initial_count, current_count));
	}
};

struct Stats {
	Stats()
		: same_color(0),
		different_color(0),
		collected(0),
		processed(0),
		collapsed(0),
		non_collapsable(0),
		cost_uncomputable(0),
		placement_uncomputable(0) {}

	std::size_t same_color;
	std::size_t different_color;
	std::size_t collected;
	std::size_t processed;
	std::size_t collapsed;
	std::size_t non_collapsable;
	std::size_t cost_uncomputable;
	std::size_t placement_uncomputable;
};

struct Visitor : CGAL::Surface_mesh_simplification::Edge_collapse_visitor_base<FaceGraph> {
	typedef GraphTraits::edges_size_type size_type;
	typedef GraphTraits::vertex_descriptor vertex_descriptor;
	typedef boost::property_map<FaceGraph, CGAL::vertex_point_t>::type Vertex_point_pmap;
	typedef CGAL::Surface_mesh_simplification::Edge_profile<FaceGraph> Profile;
	typedef boost::property_traits<Vertex_point_pmap>::value_type Point;

	Stats* stats;

	Visitor(Stats* s) : stats(s) {}

	// Called during the collecting phase for each edge collected.
	void OnCollected(Profile const&, boost::optional < Kernel::FT > const&) {
		++stats->collected;
		//std::cout << "\rEdges collected: " << stats->collected << std::flush;
	}

	// Called during the processing phase for each edge selected.
	// If cost is absent the edge won't be collapsed.
	void OnSelected(Profile const&, boost::optional<Kernel::FT> const& cost,
		size_type initial, size_type current) {
		++stats->processed;
		if (!cost)
			++stats->cost_uncomputable;

		//if (current == initial)
		//	return;
		//std::cout << "\n" << std::flush;
		//std::cout << "\rEdges remaining: " << current << std::flush;
	}

	// Called during the processing phase for each edge being collapsed.
	// If placement is absent the edge is left uncollapsed.
	void OnCollapsing(Profile const&, boost::optional<Point> const& placement) {
		if (!placement)
			++stats->placement_uncomputable;
	}

	// Called for each edge which failed the so called link-condition,
	// that is, which cannot be collapsed because doing so would
	// turn the surface mesh into a non-manifold.
	void OnNonCollapsable(Profile const&) {
		++stats->non_collapsable;
	}

	// Called AFTER each edge has been collapsed
	void OnCollapsed(Profile const&, vertex_descriptor const&) {
		++stats->collapsed;
	}
};

// BGL property map which indicates whether an edge is marked as non-removable
struct Constrained_edge_map {
	typedef boost::graph_traits<FaceGraph>::edge_descriptor key_type;

	const FaceGraph& pmesh;
	Stats* stats;
	bool use_color;
	bool use_vertex_color;
	bool has_vcolors;
	bool has_fcolors;
	int threshold;
	FaceGraph::Property_map<vertex_descriptor, CGAL::Color> vcolors;
	FaceGraph::Property_map<face_descriptor, CGAL::Color> fcolors;

	Constrained_edge_map(const FaceGraph& pmesh, Stats* stats, bool color, bool vcolor, int threshold)
		: pmesh(pmesh),
		stats(stats),
		use_color(color),
		use_vertex_color(vcolor),
		has_vcolors(false),
		has_fcolors(false),
		threshold(threshold) {
		tie(vcolors, has_vcolors) = pmesh.property_map<vertex_descriptor, CGAL::Color>("v:color");
		tie(fcolors, has_fcolors) = pmesh.property_map<face_descriptor, CGAL::Color>("f:color");
	}

	friend bool get(Constrained_edge_map map, const key_type& edge) {
		if (!map.use_color)
			return false;

		auto v1 = map.pmesh.vertex(edge, 0);
		auto v2 = map.pmesh.vertex(edge, 1);
		CGAL::Color c1;
		CGAL::Color c2;

		if (map.has_fcolors && !map.use_vertex_color) {
			auto f1 = map.pmesh.face(map.pmesh.halfedge(v1));
			auto f2 = map.pmesh.face(map.pmesh.halfedge(v2));
			c1 = map.fcolors[f1];
			c2 = map.fcolors[f2];
		} else if (map.has_vcolors) {
			c1 = map.vcolors[v1];
			c2 = map.vcolors[v2];
		}

		if (c1 == c2)
			map.stats->same_color++;
		else
			map.stats->different_color++;

		int r = c1.red() - c2.red();
		int g = c1.green() - c2.green();
		int b = c1.blue() - c2.blue();
		return r * r + g * g + b * b > map.threshold*map.threshold;
	}
};

class ColorConstrainedSimplification {
public:
	typedef FaceGraph TM;

	typedef boost::graph_traits<TM> GraphTraits;
	typedef GraphTraits::vertex_descriptor      vertex_descriptor;
	typedef GraphTraits::vertex_iterator        vertex_iterator;
	typedef GraphTraits::halfedge_descriptor    halfedge_descriptor;
	typedef GraphTraits::halfedge_iterator      halfedge_iterator;
	typedef CGAL::Halfedge_around_source_iterator<TM> out_edge_iterator;
	typedef CGAL::Halfedge_around_target_iterator<TM> in_edge_iterator;
	typedef GraphTraits::traversal_category     traversal_category;
	typedef GraphTraits::edges_size_type        size_type;
	typedef GraphTraits::edge_iterator          edge_iterator;

	
	typedef boost::lazy_disable_if<boost::is_const<CGAL::Point_3<CGAL::Epick>>, CGAL::internal::
		Get_vertex_point_map_for_Surface_mesh_return_type<CGAL::Point_3<CGAL::Epick>>>::type VertexPointMap;
	typedef CGAL::Surface_mesh_simplification::Edge_profile<TM, VertexPointMap> Profile;
	typedef VertexPointMap Vertex_point_pmap;
	typedef boost::property_traits<Vertex_point_pmap>::value_type Point;

	typedef Kernel Traits;
	typedef Traits::Equal_3 Equal_3;
	typedef Traits::Vector_3 Vector;
	typedef Traits::FT FT;

	typedef boost::optional<FT> Cost_type;
	typedef boost::optional<Point> Placement_type;
	typedef ColorConstrainedSimplification Self;

	typedef Constrained_edge_map EdgeIsConstrainedMap;

	struct Compare_id {
		Compare_id() : mAlgorithm(0) {}

		Compare_id(Self const* aAlgorithm) : mAlgorithm(aAlgorithm) {}

		bool operator() (halfedge_descriptor const& a, halfedge_descriptor const& b) const {
			return mAlgorithm->get_halfedge_id(a) < mAlgorithm->get_halfedge_id(b);
		}

		Self const* mAlgorithm;
	};

	struct Compare_cost {
		Compare_cost() : mAlgorithm(nullptr) {}

		Compare_cost(Self const* aAlgorithm) : mAlgorithm(aAlgorithm) {}

		bool operator() (halfedge_descriptor const& a, halfedge_descriptor const& b) const {
			// NOTE: A cost is an optional<> value.
			// Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
			// In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
			return mAlgorithm->get_data(a).cost() < mAlgorithm->get_data(b).cost();
		}

		Self const* mAlgorithm;
	};

	struct edge_id : boost::put_get_helper<size_type, edge_id> {
		edge_id() : mAlgorithm(nullptr) {}

		edge_id(Self const* aAlgorithm) : mAlgorithm(aAlgorithm) {}

		size_type operator[] (halfedge_descriptor const& e) const {
			return mAlgorithm->get_edge_id(e);
		}

		Self const* mAlgorithm;
	};

	typedef CGAL::Modifiable_priority_queue<halfedge_descriptor, Compare_cost, edge_id> PQ;
	typedef PQ::handle pq_handle;

	// An Edge_data is associated with EVERY _ edge in the mesh (collapseable or not).
	// It relates the edge with the PQ-handle needed to update the priority queue
	// It also relates the edge with a policy-based cache
	class Edge_data {
	public:
		Edge_data() : mPQHandle() {}

		Cost_type const& cost() const { return mCost; }
		Cost_type      & cost() { return mCost; }

		pq_handle PQ_handle() const { return mPQHandle; }

		bool is_in_PQ() const { return mPQHandle != PQ::null_handle(); }

		void set_PQ_handle(pq_handle h) { mPQHandle = h; }

		void reset_PQ_handle() { mPQHandle = PQ::null_handle(); }

	private:
		Cost_type mCost;
		pq_handle mPQHandle;
	};

	typedef Edge_data* Edge_data_ptr;
	typedef boost::scoped_array<Edge_data> Edge_data_array;

	ColorConstrainedSimplification(TM& aSurface
		, Custom_stop_predicate& aShould_stop
//		, VertexIndexMap        const& aVertex_index_map
//		, VertexPointMap        const& aVertex_point_map
//		, EdgeIndexMap          const& aEdge_index_map
		, EdgeIsConstrainedMap& aEdge_is_constrained_map
//		, CGAL::Surface_mesh_simplification::LindstromTurk_cost<TM> const& aGet_cost
//		, CGAL::Surface_mesh_simplification::LindstromTurk_placement<TM> const& aGet_placement
//		, VisitorT              aVisitor
	)
		: mSurface(aSurface)
		, Should_stop(aShould_stop)
		, Vertex_index_map(get(boost::vertex_index, mSurface))
		, Vertex_point_map(get(boost::vertex_point, mSurface))
		, Edge_index_map(get(boost::halfedge_index, mSurface))
		, Edge_is_constrained_map(aEdge_is_constrained_map)
		, Get_cost(CGAL::Surface_mesh_simplification::LindstromTurk_cost<TM>())
		, Get_placement(CGAL::Surface_mesh_simplification::LindstromTurk_placement<TM>())
//		, m_has_border(false)
	{
		const FT cMaxDihedralAngleCos = std::cos(1.0 * CGAL_PI / 180.0);
		mcMaxDihedralAngleCos2 = cMaxDihedralAngleCos * cMaxDihedralAngleCos;
		//Get_cost = CGAL::Surface_mesh_simplification::LindstromTurk_cost<TM>();
		//Get_placement = CGAL::Surface_mesh_simplification::LindstromTurk_placement<TM>();
		//Vertex_point_map = get(boost::vertex_point, mSurface);
		//Edge_index_map = get(boost::halfedge_index, mSurface);

		halfedge_iterator eb, ee;
		for (boost::tie(eb, ee) = halfedges(mSurface); eb != ee; ++eb) {
			halfedge_descriptor ed = *eb;
			if (is_border(ed)) {
				m_has_border = true;
				break;
			}
		}
	}

	void simplify() {
		// Loop over all the _undirected_ edges in the surface putting them in the PQ
		Equal_3 equal_points = Traits().equal_3_object();
		size_type lSize = num_edges(mSurface);
		mInitialEdgeCount = mCurrentEdgeCount = static_cast<size_type>(
			std::distance(boost::begin(edges(mSurface)), boost::end(edges(mSurface))));;

		mEdgeDataArray.reset(new Edge_data[lSize]);

		mPQ.reset(new PQ(lSize, Compare_cost(this), edge_id(this)));

		std::size_t id = 0;
		std::set<halfedge_descriptor> zero_length_edges;

		edge_iterator eb, ee;
		for (boost::tie(eb, ee) = edges(mSurface); eb != ee; ++eb, id += 2) {
			halfedge_descriptor lEdge = halfedge(*eb, mSurface);

			if (is_constrained(lEdge))
				continue; //no not insert constrained edges

			Profile const& lProfile = create_profile(lEdge);
			if (!equal_points(lProfile.p0(), lProfile.p1())) {
				Edge_data& lData = get_data(lEdge);

				lData.cost() = get_cost(lProfile);
				insert_in_PQ(lEdge, lData);
			} else {
				zero_length_edges.insert(primary_edge(lEdge));
			}
		}

		for (auto it = zero_length_edges.begin(), it_end = zero_length_edges.end(); it != it_end; ++it) {
			Profile const& lProfile = create_profile(*it);

			if (!Is_collapse_topologically_valid(lProfile))
				continue;

			// edges of length 0 removed no longer need to be treated
			if (lProfile.left_face_exists()) {
				halfedge_descriptor lEdge_to_remove = is_constrained(lProfile.vL_v0()) ?
					primary_edge(lProfile.v1_vL()) :
					primary_edge(lProfile.vL_v0());
				zero_length_edges.erase(lEdge_to_remove);
				Edge_data& lData = get_data(lEdge_to_remove);
				if (lData.is_in_PQ()) {
					remove_from_PQ(lEdge_to_remove, lData);
				}
				--mCurrentEdgeCount;
			}

			if (lProfile.right_face_exists()) {
				halfedge_descriptor lEdge_to_remove = is_constrained(lProfile.vR_v1()) ?
					primary_edge(lProfile.v0_vR()) :
					primary_edge(lProfile.vR_v1());
				zero_length_edges.erase(lEdge_to_remove);
				Edge_data& lData = get_data(lEdge_to_remove);
				if (lData.is_in_PQ()) {
					remove_from_PQ(lEdge_to_remove, lData);
				}
				--mCurrentEdgeCount;
			}

			--mCurrentEdgeCount;

			//the placement is trivial, it's always the point itself
			Placement_type lPlacement = lProfile.p0();
			vertex_descriptor rResult = halfedge_collapse_bk_compatibility(lProfile.v0_v1(), Edge_is_constrained_map);
			boost::put(Vertex_point_map, rResult, *lPlacement);
		}

		// Pops and processes each edge from the PQ
		boost::optional<halfedge_descriptor> lEdge;
		while ((lEdge = pop_from_PQ())) {
			Profile const& lProfile = create_profile(*lEdge);
			Cost_type lCost = get_data(*lEdge).cost();

			if (lCost) {
				if (Should_stop(*lCost, lProfile, mInitialEdgeCount, mCurrentEdgeCount)) {
					//Visitor.OnStopConditionReached(lProfile);
					break;
				}

				if (Is_collapse_topologically_valid(lProfile)) {
					// The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
					// satisfying its constraints. In that case the remaining vertex is simply left unmoved.
					Placement_type lPlacement = get_placement(lProfile);

					if (Is_collapse_geometrically_valid(lProfile, lPlacement)) {
						Collapse(lProfile, lPlacement);
					}
				} else {
					//Visitor.OnNonCollapsable(lProfile);
				}
			}
		}
	}
private:
	void Collapse(Profile const& aProfile, Placement_type aPlacement) {
		//Visitor.OnCollapsing(aProfile, aPlacement);
		--mCurrentEdgeCount;

		// If the top/bottom facets exists, they are removed and the edges v0vt and Q-B along with them.
		// In that case their corresponding pairs must be pop off the queue
		if (aProfile.left_face_exists()) {
			halfedge_descriptor lV0VL = primary_edge(aProfile.vL_v0());
			if (is_constrained(lV0VL)) //make sure a constrained edge will not disappear
				lV0VL = primary_edge(aProfile.v1_vL());

			Edge_data& lData = get_data(lV0VL);
			if (lData.is_in_PQ()) {
				remove_from_PQ(lV0VL, lData);
			}
			--mCurrentEdgeCount;
		}
		if (aProfile.right_face_exists()) {
			halfedge_descriptor lVRV1 = primary_edge(aProfile.vR_v1());
			if (is_constrained(lVRV1)) //make sure a constrained edge will not disappear
				lVRV1 = primary_edge(aProfile.v0_vR());

			Edge_data& lData = get_data(lVRV1);
			if (lData.is_in_PQ()) {
				remove_from_PQ(lVRV1, lData);
			}
			--mCurrentEdgeCount;
		}
		// Perform the actual collapse.
		// This is an external function.
		// It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ,PT and QB
		// (PT and QB are removed if they are not null).
		// All other edges must be kept.
		// All directed edges incident to vertex removed are relink to the vertex kept.
		vertex_descriptor rResult = halfedge_collapse_bk_compatibility(aProfile.v0_v1(), Edge_is_constrained_map);

		if (aPlacement) {
			boost::put(Vertex_point_map, rResult, *aPlacement);
		}
		//Visitor.OnCollapsed(aProfile, rResult);
		Update_neighbors(rResult);
	}

	void Update_neighbors(vertex_descriptor const& aKeptV) {
		// (A) Collect all edges to update their cost: all those around each vertex adjacent to the vertex kept
		typedef std::set<halfedge_descriptor, Compare_id> edges;

		edges lToUpdate(Compare_id(this));
		edges lToInsert(Compare_id(this));

		// (A.1) Loop around all vertices adjacent to the vertex kept
		in_edge_iterator eb1, ee1;
		for (boost::tie(eb1, ee1) = halfedges_around_target(aKeptV, mSurface); eb1 != ee1; ++eb1) {
			halfedge_descriptor lEdge1 = *eb1;
			vertex_descriptor lAdj_k = source(lEdge1, mSurface);

			// (A.2) Loop around all edges incident on each adjacent vertex
			in_edge_iterator eb2, ee2;
			for (boost::tie(eb2, ee2) = halfedges_around_target(lAdj_k, mSurface); eb2 != ee2; ++eb2) {
				halfedge_descriptor lEdge2 = primary_edge(*eb2);
				Edge_data& lData2 = get_data(lEdge2);

				// Only edges still in the PQ needs to be updated, the other needs to be re-inserted
				if (lData2.is_in_PQ())
					lToUpdate.insert(lEdge2);
				else
					lToInsert.insert(lEdge2);
			}
		}
		// (B) Proceed to update the costs.
		for (auto it = lToUpdate.begin(), eit = lToUpdate.end(); it != eit; ++it) {
			halfedge_descriptor lEdge = *it;
			Edge_data& lData = get_data(lEdge);
			Profile const& lProfile = create_profile(lEdge);
			lData.cost() = get_cost(lProfile);
			update_in_PQ(lEdge, lData);
		}

		// (C) Insert ignored edges
		// I think that this should be done for edges eliminated because of the geometric criteria
		// and not the topological one.However maintaining such a set might be more expensive
		// and hard to be safe ...
		for (auto it = lToInsert.begin(), eit = lToInsert.end(); it != eit; ++it) {
			halfedge_descriptor lEdge = *it;
			if (is_constrained(lEdge)) continue; //do not insert constrained edges
			Edge_data& lData = get_data(lEdge);
			Profile const& lProfile = create_profile(lEdge);
			lData.cost() = get_cost(lProfile);
			insert_in_PQ(lEdge, lData);
		}
	}

	bool is_primary_edge(halfedge_descriptor const& aEdge) const {
		return (get_halfedge_id(aEdge) % 2) == 0;
	}

	Cost_type get_cost(Profile const& aProfile) const {
		return Get_cost(aProfile, get_placement(aProfile));
	}

	Placement_type get_placement(Profile const& aProfile) const {
		return Get_placement(aProfile);
	}

	halfedge_descriptor primary_edge(halfedge_descriptor const& aEdge) {
		return is_primary_edge(aEdge) ? aEdge : opposite(aEdge, mSurface);
	}

	Profile create_profile(halfedge_descriptor const& aEdge) {
		return Profile(aEdge, mSurface, Vertex_index_map, Vertex_point_map, Edge_index_map, m_has_border);
	}

	bool is_constrained(halfedge_descriptor const& aEdge) const {
		return get(Edge_is_constrained_map, edge(aEdge, mSurface));
	}

	size_type get_halfedge_id(halfedge_descriptor const& aEdge) const {
		return Edge_index_map[aEdge];
	}

	size_type get_edge_id(halfedge_descriptor const& aEdge) const {
		return get_halfedge_id(aEdge) / 2;
	}

	Edge_data& get_data(halfedge_descriptor const& aEdge) const {
		return mEdgeDataArray[get_edge_id(aEdge)];
	}

	void insert_in_PQ(halfedge_descriptor const& aEdge, Edge_data& aData) {
		aData.set_PQ_handle(mPQ->push(aEdge));
	}

	void update_in_PQ(halfedge_descriptor const& aEdge, Edge_data& aData) {
		aData.set_PQ_handle(mPQ->update(aEdge, aData.PQ_handle()));
	}

	void remove_from_PQ(halfedge_descriptor const& aEdge, Edge_data& aData) {
		aData.set_PQ_handle(mPQ->erase(aEdge, aData.PQ_handle()));
	}

	boost::optional<halfedge_descriptor> pop_from_PQ() {
		boost::optional<halfedge_descriptor> rEdge = mPQ->extract_top();
		if (rEdge) {
			get_data(*rEdge).reset_PQ_handle();
		}
		return rEdge;
	}

	bool is_border(vertex_descriptor const& aV) const {
		bool rR = false;

		in_edge_iterator eb, ee;
		for (boost::tie(eb, ee) = halfedges_around_target(aV, mSurface); eb != ee; ++eb) {
			halfedge_descriptor lEdge = *eb;
			if (is_border(lEdge)) {
				rR = true;
				break;
			}
		}

		return rR;
	}

	bool is_edge_a_border(halfedge_descriptor const& aEdge) const {
		return is_border(aEdge) || is_border(opposite(aEdge, mSurface));
	}

	bool is_border(halfedge_descriptor const& aEdge) const {
		return face(aEdge, mSurface) == boost::graph_traits<TM>::null_face();
	}

	bool is_constrained(vertex_descriptor const& aV) const {
		in_edge_iterator eb, ee;
		for (boost::tie(eb, ee) = halfedges_around_target(aV, mSurface); eb != ee; ++eb)
			if (is_constrained(*eb)) return true;
		return false;
	}

	// Some edges are NOT collapse able: doing so would break the topological consistency of the mesh.
	// This function returns true if a edge 'p->q' can be collapsed.
	//
	// An edge p->q can be collapsed iff it satisfies the "link condition"
	// (as described in the "Mesh Optimization" article of Hoppe et al (1993))
	//
	// The link condition is as follows: for every vertex 'k' adjacent to both 'p and 'q',
	// "p,k,q" is a facet of the mesh.
	//
	bool Is_collapse_topologically_valid(Profile const& aProfile) {
		bool rR = true;
		out_edge_iterator eb1, ee1;
		out_edge_iterator eb2, ee2;

		// Simple tests handling the case of non-manifold situations at a vertex or edge (pinching)
		// (even if we advertise one should not use a surface mesh with such features)
		if (aProfile.left_face_exists()) {
			if (CGAL::is_border(opposite(aProfile.v1_vL(), mSurface), mSurface) &&
				CGAL::is_border(opposite(aProfile.vL_v0(), mSurface), mSurface)
				) return false;

			if (aProfile.right_face_exists() &&
				CGAL::is_border(opposite(aProfile.vR_v1(), mSurface), mSurface) &&
				CGAL::is_border(opposite(aProfile.v0_vR(), mSurface), mSurface)
				) return false;
		} else {
			if (aProfile.right_face_exists()) {
				if (CGAL::is_border(opposite(aProfile.vR_v1(), mSurface), mSurface) &&
					CGAL::is_border(opposite(aProfile.v0_vR(), mSurface), mSurface)
					) return false;
			} else
				return false;
		}

		// The following loop checks the link condition for v0_v1.
		// Specifically, that for every vertex 'k' adjacent to both 'p and 'q', 'pkq' is a face of the mesh.
		for (boost::tie(eb1, ee1) = halfedges_around_source(aProfile.v0(), mSurface); rR && eb1 != ee1; ++eb1) {
			halfedge_descriptor v0_k = *eb1;

			if (v0_k != aProfile.v0_v1()) {
				vertex_descriptor k = target(v0_k, mSurface);

				for (boost::tie(eb2, ee2) = halfedges_around_source(k, mSurface); rR && eb2 != ee2; ++eb2) {
					halfedge_descriptor k_v1 = *eb2;

					if (target(k_v1, mSurface) == aProfile.v1()) {
						// At this point we know p-q-k are connected and we need to determine if this triangle is a face of the mesh.
						//
						// Since the mesh is known to be triangular there are at most two faces sharing the edge p-q.
						//
						// If p->q is NOT a border edge, the top face is p->q->t where t is target(next(p->q))
						// If q->p is NOT a border edge, the bottom face is q->p->b where b is target(next(q->p))
						//
						// If k is either t or b then p-q-k *might* be a face of the mesh. It won't be if k==t but p->q is border
						// or k==b but q->b is a border (because in that case even though there exists triangles p->q->t (or q->p->b)
						// they are holes, not faces)
						bool lIsFace = (aProfile.vL() == k && aProfile.left_face_exists())
							|| (aProfile.vR() == k && aProfile.right_face_exists());

						if (!lIsFace) {
							rR = false;
							break;
						}
					}
				}
			}
		}

		if (rR) {
			// ensure two constrained edges cannot get merged
			if (is_edge_adjacent_to_a_constrained_edge(aProfile, Edge_is_constrained_map))
				return false;

			if (aProfile.is_v0_v1_a_border()) {
				if (Is_open_triangle(aProfile.v0_v1())) {
					rR = false;
				}
			} else if (aProfile.is_v1_v0_a_border()) {
				if (Is_open_triangle(aProfile.v1_v0())) {
					rR = false;
				}
			} else {
				if (is_border(aProfile.v0()) && is_border(aProfile.v1())) {
					rR = false;
				} else {
					bool lTetra = is_tetrahedron(aProfile.v0_v1(), mSurface);

					if (lTetra) {
						rR = false;
					}

					if (next(aProfile.v0_v1(), mSurface) == opposite(prev(aProfile.v1_v0(), mSurface), mSurface) &&
						prev(aProfile.v0_v1(), mSurface) == opposite(next(aProfile.v1_v0(), mSurface), mSurface)) {
						return false;
					}
				}
			}
		}
		return rR;
	}

	halfedge_descriptor find_connection(vertex_descriptor const& v0, vertex_descriptor const& v1) const {
		out_edge_iterator eb, ee;
		for (boost::tie(eb, ee) = halfedges_around_source(v0, mSurface); eb != ee; ++eb) {
			halfedge_descriptor out = *eb;
			if (target(out, mSurface) == v1)
				return out;
		}
		return halfedge_descriptor();
	}

	bool Is_collapse_geometrically_valid(Profile const& aProfile, Placement_type k0) {
		bool rR = true;
		if (k0) {
			// Use the current link to extract all local triangles incident to 'vx' in the collapsed mesh (which at this point doesn't exist yet)
			typedef Profile::vertex_descriptor_vector::const_iterator link_iterator;
			link_iterator linkb = aProfile.link().begin();
			link_iterator linke = aProfile.link().end();
			link_iterator linkl = prev(linke);

			for (link_iterator l = linkb; l != linke && rR; ++l) {
				link_iterator pv = (l == linkb ? linkl : prev(l));
				link_iterator nx = (l == linkl ? linkb : next(l));

				// k0,k1 and k3 are three consecutive vertices along the link.
				vertex_descriptor k1 = *pv;
				vertex_descriptor k2 = *l;
				vertex_descriptor k3 = *nx;
				halfedge_descriptor e12 = find_connection(k1, k2);
				halfedge_descriptor e23 = k3 != k1 ? find_connection(k2, k3) : halfedge_descriptor();

				// If 'k1-k2-k3' are connected there will be two adjacent triangles 'k0,k1,k2' and 'k0,k2,k3' after the collapse.
				if (handle_assigned(e12) && handle_assigned(e23)) {
					if (!are_shared_triangles_valid(*k0, get_point(k1), get_point(k2), get_point(k3))) {
						rR = false;
					}
				}

				if (rR) {
					// Also check the triangles 'k0,k1,k2' and it's adjacent along e12: 'k4,k2,k1', if exist
					vertex_descriptor k4 = find_exterior_link_triangle_3rd_vertex(e12, aProfile.v0(), aProfile.v1());

					// There is indeed a triangle shared along e12
					if (handle_assigned(k4)) {
						if (!are_shared_triangles_valid(get_point(k1), get_point(k4), get_point(k2), *k0)) {
							rR = false;
						}
					}
				}

				if (rR) {
					// And finally, check the triangles 'k0,k2,k3' and it's adjacent e23: 'k5,k3,k2' if exist
					vertex_descriptor k5 = find_exterior_link_triangle_3rd_vertex(e23, aProfile.v0(), aProfile.v1());

					// There is indeed a triangle shared along e12
					if (handle_assigned(k5)) {
						if (!are_shared_triangles_valid(get_point(k2), get_point(k5), get_point(k3), *k0)) {
							rR = false;
						}
					}
				}
			}
		}
		return rR;
	}

	bool are_shared_triangles_valid(Point const& p0, Point const& p1, Point const& p2, Point const& p3) const {
		bool rR = false;

		Vector e01 = Traits().construct_vector_3_object()(p0, p1);
		Vector e02 = Traits().construct_vector_3_object()(p0, p2);
		Vector e03 = Traits().construct_vector_3_object()(p0, p3);

		Vector n012 = Traits().construct_cross_product_vector_3_object()(e01, e02);
		Vector n023 = Traits().construct_cross_product_vector_3_object()(e02, e03);

		FT l012 = Traits().compute_scalar_product_3_object()(n012, n012);
		FT l023 = Traits().compute_scalar_product_3_object()(n023, n023);

		FT larger = (std::max)(l012, l023);
		FT smaller = (std::min)(l012, l023);

		const double cMaxAreaRatio = 1e8;

		if (larger < cMaxAreaRatio * smaller) {
			FT l0123 = Traits().compute_scalar_product_3_object()(n012, n023);

			if (CGAL_NTS is_positive(l0123)) {
				rR = true;
			} else {
				if ((l0123 * l0123) <= mcMaxDihedralAngleCos2 * (l012 * l023)) {
					rR = true;
				}
			}
		}

		return rR;
	}

	vertex_descriptor find_exterior_link_triangle_3rd_vertex(halfedge_descriptor const& e, vertex_descriptor const& v0, vertex_descriptor const& v1) const {
		vertex_descriptor r;
		if (handle_assigned(e)) {
			vertex_descriptor ra = target(next(e, mSurface), mSurface);
			vertex_descriptor rb = source(prev(e, mSurface), mSurface);

			if (ra == rb && ra != v0 && ra != v1) {
				r = ra;
			} else {
				ra = target(next(opposite(e, mSurface), mSurface), mSurface);
				rb = source(prev(opposite(e, mSurface), mSurface), mSurface);

				if (ra == rb && ra != v0 && ra != v1) {
					r = ra;
				}
			}
		}
		return r;
	}

	template<class Handle>
	bool handle_assigned(Handle h) const {
		Handle null;
		return h != null;
	}

	boost::property_traits<VertexPointMap>::reference get_point(vertex_descriptor const& aV) const {
		return get(Vertex_point_map, aV);
	}

	bool Is_open_triangle(halfedge_descriptor const& h1) {
		bool rR = false;
		halfedge_descriptor h2 = next(h1, mSurface);
		halfedge_descriptor h3 = next(h2, mSurface);

		// First check if it is a triangle 
		if (next(h3, mSurface) == h1) {
			rR = is_border(h2) && is_border(h3);
		}

		return rR;
	}

	vertex_descriptor halfedge_collapse_bk_compatibility(halfedge_descriptor const& pq, EdgeIsConstrainedMap aEdge_is_constrained_map) {
		return CGAL::Euler::collapse_edge(edge(pq, mSurface), mSurface, aEdge_is_constrained_map);
	}

	bool is_edge_adjacent_to_a_constrained_edge(Profile const& aProfile, EdgeIsConstrainedMap) {
		return is_constrained(aProfile.v0()) && is_constrained(aProfile.v1());
	}

	TM& mSurface;
	Custom_stop_predicate& Should_stop;
	CGAL::SM_index_pmap<CGAL::Point_3<CGAL::Epick>, CGAL::SM_Vertex_index>& Vertex_index_map;
	VertexPointMap& Vertex_point_map;
	CGAL::SM_index_pmap<CGAL::Point_3<CGAL::Epick>, CGAL::SM_Halfedge_index>& Edge_index_map;
	EdgeIsConstrainedMap& Edge_is_constrained_map;
	CGAL::Surface_mesh_simplification::LindstromTurk_cost<TM>& Get_cost;
	CGAL::Surface_mesh_simplification::LindstromTurk_placement<TM>& Get_placement;
	//VisitorT Visitor;
	bool m_has_border;
	Edge_data_array mEdgeDataArray;
	boost::scoped_ptr<PQ> mPQ;
	size_type mInitialEdgeCount;
	size_type mCurrentEdgeCount;
	FT mcMaxDihedralAngleCos2;
};

class Polyhedron_demo_mesh_simplification_color_plugin :
	public QObject,
	public CGAL::Three::Polyhedron_demo_plugin_interface {
	Q_OBJECT
		Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
		Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
	QList<QAction*> actions() const override {
		return _actions;
	}

	void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) override {
		mw = mainWindow;
		scene = scene_interface;
		QAction *actionSimplify = new QAction("Color Constrained Simplification", mw);
		actionSimplify->setProperty("subMenuName", "Triangulated Surface Mesh Simplification");
		connect(actionSimplify, SIGNAL(triggered()), this, SLOT(on_actionSimplify_triggered()));
		_actions << actionSimplify;
	}

	bool applicable(QAction*) const override {
		return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()))
			|| qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
	}
public Q_SLOTS:
	void on_actionSimplify_triggered();
private:
	CGAL::Three::Scene_interface *scene;
	QMainWindow *mw;
	QList<QAction*> _actions;

}; // end Polyhedron_demo_mesh_simplification_color_plugin

void Polyhedron_demo_mesh_simplification_color_plugin::on_actionSimplify_triggered() {
	const int index = scene->mainSelectionIndex();
	auto poly_item = qobject_cast<Scene_facegraph_item*>(scene->item(index));
	auto selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

	if (poly_item || selection_item) {
		FaceGraph& pmesh = poly_item != nullptr ? *poly_item->polyhedron() : *selection_item->polyhedron();

		// get option
		QDialog dialog(mw);
		Ui::Mesh_simplification_color_dialog ui;
		ui.setupUi(&dialog);
		connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
		connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

		auto bbox = poly_item != nullptr ? poly_item->bbox()
			: selection_item != nullptr ? selection_item->bbox() : scene->bbox();

		double diago_length = CGAL::sqrt((bbox.xmax() - bbox.xmin())*(bbox.xmax() - bbox.xmin())
			+ (bbox.ymax() - bbox.ymin())*(bbox.ymax() - bbox.ymin()) +
			(bbox.zmax() - bbox.zmin())*(bbox.zmax() - bbox.zmin()));

		ui.m_nb_edges->setValue((int)(num_halfedges(pmesh) / 4));
		ui.m_nb_edges->setMaximum((int)num_halfedges(pmesh));
		ui.m_edge_length->setValue(diago_length * 0.05);
		ui.m_color_threshold->setValue(50);

		bool has_vcolors;
		bool has_fcolors;
		FaceGraph::Property_map<vertex_descriptor, CGAL::Color> vcolors;
		FaceGraph::Property_map<face_descriptor, CGAL::Color> fcolors;
		tie(vcolors, has_vcolors) = pmesh.property_map<vertex_descriptor, CGAL::Color>("v:color");
		tie(fcolors, has_fcolors) = pmesh.property_map<face_descriptor, CGAL::Color>("f:color");

		if (!has_fcolors && has_vcolors) {
			std::cout << "\nAdding face colors..." << std::endl;
			pmesh.add_property_map<face_descriptor, CGAL::Color>("f:color", CGAL::Color(0, 0, 0));
			tie(fcolors, has_fcolors) = pmesh.property_map<face_descriptor, CGAL::Color>("f:color");

			BOOST_FOREACH(FaceGraph::Face_index f, pmesh.faces()) {
				unsigned int r = 0, g = 0, b = 0;
				BOOST_FOREACH(FaceGraph::Vertex_index v, CGAL::vertices_around_face(pmesh.halfedge(f), pmesh)) {
					CGAL::Color c = vcolors[v];
					r += c.red() * c.red();
					g += c.green() * c.green();
					b += c.blue() * c.blue();
				}
				CGAL::Color fcolor(
					static_cast<unsigned char>(sqrt(r / 3)),
					static_cast<unsigned char>(sqrt(g / 3)),
					static_cast<unsigned char>(sqrt(b / 3))
				);
				fcolors[f] = fcolor;
			}
			std::cout << "Face colors were added..." << std::endl;
		}

		if (has_vcolors)
			ui.m_source->addItem("Vertex");
		if (has_fcolors)
			ui.m_source->addItem("Face");
		if (!has_vcolors && !has_fcolors) {
			ui.m_source->setEnabled(false);
			ui.m_use_source->setChecked(false);
			ui.m_use_source->setEnabled(false);
		}

		std::cout << "\nHas vertex color: " << std::boolalpha << has_vcolors << std::endl
			<< "Has face color: " << has_fcolors << std::endl;

		// check user cancellation
		if (dialog.exec() == QDialog::Rejected)
			return;

		// simplify
		QTime time;
		time.start();
		QApplication::setOverrideCursor(Qt::WaitCursor);
		QApplication::processEvents();

		Stats stats;
		Visitor visitor(&stats);
		Constrained_edge_map bem(pmesh, &stats,
			ui.m_use_source->isChecked(),
			ui.m_source->currentIndex() == 0,
			ui.m_color_threshold->value()
		);
		Custom_stop_predicate stop(
			ui.m_combinatorial->currentIndex() == 0 && !ui.m_use_nb_edges->isChecked() && !ui.m_use_edge_length->isChecked(),
			ui.m_use_nb_edges->isChecked() ? ui.m_nb_edges->value() : 0,
			ui.m_use_edge_length->isChecked() ? ui.m_edge_length->value() : std::numeric_limits<double>::max()
		);

		if (ui.m_base_placement->currentIndex() == 0) {
			typedef SMS::LindstromTurk_placement<FaceGraph> BasePlacement;

			if (ui.m_use_bounded_normal_change_placement->isChecked()) {
				typedef SMS::Bounded_normal_change_placement<BasePlacement> Placement;
				typedef SMS::Constrained_placement<Placement, Constrained_edge_map> ConstrainedPlacement;
				edge_collapse(pmesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
					.get_placement(ConstrainedPlacement(bem))
					.visitor(visitor));
			} else {
				typedef BasePlacement Placement;
				typedef SMS::Constrained_placement<Placement, Constrained_edge_map> ConstrainedPlacement;
				edge_collapse(pmesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
					.get_placement(ConstrainedPlacement(bem))
					.visitor(visitor));
			}
		} else if (ui.m_base_placement->currentIndex() == 1) {
			typedef SMS::Midpoint_placement<FaceGraph> BasePlacement;

			if (ui.m_use_bounded_normal_change_placement->isChecked()) {
				typedef SMS::Bounded_normal_change_placement<BasePlacement> Placement;
				typedef SMS::Constrained_placement<Placement, Constrained_edge_map> ConstrainedPlacement;
				edge_collapse(pmesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
					.get_placement(ConstrainedPlacement(bem))
					.visitor(visitor));
			} else {
				typedef BasePlacement Placement;
				typedef SMS::Constrained_placement<Placement, Constrained_edge_map> ConstrainedPlacement;
				edge_collapse(pmesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
					.get_placement(ConstrainedPlacement(bem))
					.visitor(visitor));
			}
		}

		if (ui.m_use_source->isChecked()) {
			std::cout << "\nSame color: " << stats.same_color << std::endl
				<< "Different color: " << stats.different_color << std::endl;
		}
		std::cout << "\nEdges collected: " << stats.collected << std::endl
			<< "Edges processed: " << stats.processed << std::endl
			<< "Edges collapsed: " << stats.collapsed << std::endl
			<< "Edges not collapsed due to topological constraints: " << stats.non_collapsable << std::endl
			<< "Edge not collapsed due to cost computation constraints: " << stats.cost_uncomputable << std::endl
			<< "Edge not collapsed due to placement computation constraints: " << stats.placement_uncomputable << std::endl
			<< "Time elapsed: " << time.elapsed() << " ms" << std::endl;

		// update scene
		if (poly_item != nullptr) {
			poly_item->invalidateOpenGLBuffers();
			poly_item->polyhedron()->collect_garbage();
		} else {
			selection_item->polyhedron_item()->polyhedron()->collect_garbage();
			selection_item->poly_item_changed();
			selection_item->changed_with_poly_item();
			selection_item->invalidateOpenGLBuffers();
		}
		scene->itemChanged(index);
		QApplication::restoreOverrideCursor();
	}
}
#include "Mesh_simplification_color_plugin.moc"