#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "ui_Mesh_simplification_color_dialog.h"

#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>

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
		placement_uncomputable(0)
	{}

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
		}
		else if (map.has_vcolors) {
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
		return r*r + g*g + b*b > map.threshold*map.threshold;
	}
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

		ui.m_nb_edges->setValue((int) (num_halfedges(pmesh) / 4));
		ui.m_nb_edges->setMaximum((int) num_halfedges(pmesh));
		ui.m_edge_length->setValue(diago_length * 0.05);
		ui.m_color_threshold->setValue(50);
		
		bool has_vcolors = pmesh.property_map<vertex_descriptor, CGAL::Color>("v:color").second;
		bool has_fcolors = pmesh.property_map<face_descriptor, CGAL::Color>("f:color").second;
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
			}
			else {
				typedef BasePlacement Placement;
				typedef SMS::Constrained_placement<Placement, Constrained_edge_map> ConstrainedPlacement;
				edge_collapse(pmesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
					.get_placement(ConstrainedPlacement(bem))
					.visitor(visitor));
			}
		}
		else if (ui.m_base_placement->currentIndex() == 1) {
			typedef SMS::Midpoint_placement<FaceGraph> BasePlacement;

			if (ui.m_use_bounded_normal_change_placement->isChecked()) {
				typedef SMS::Bounded_normal_change_placement<BasePlacement> Placement;
				typedef SMS::Constrained_placement<Placement, Constrained_edge_map> ConstrainedPlacement;
				edge_collapse(pmesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
					.get_placement(ConstrainedPlacement(bem))
					.visitor(visitor));
			}
			else {
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
		}
		else {
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