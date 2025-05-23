#include "facet_classification.h"
#include "meshIO.h"
#include "CLI.hpp"
#include "meshAlgorithm.h"
#include "remesh.h"
#include "SurfaceHoleFilling.h"
#include <iostream>
#include "igl/bfs_orient.h"

#define _DEBUG_ 1

using namespace std;
using std::cout;

int main(int argc, char** argv)
{
	CLI::App app{ "MeshConveter" };
	std::cout << "Program version: 2023.11.22" << std::endl;
	vector<string> input_filenames;
	string input_filename_ex;
	bool exportMESH = false;
	bool exportVTK = false;
	bool exportEpsVTK = false;
	bool exportPLY = false;
	bool exportPLS = false;
	bool exportFacet = false;
	bool exportOBJ = false;
	bool exportStlIn = false;
	bool exportSTL = false;
	bool shuffleMark = false;
	bool resetOritation = false;
	bool checkOritation = false;
	bool reverseFacetOrient = false;
	bool RepairZeroAera = false;
	bool resetOritationFaceid = false;
	bool removebox = false;
	bool bremesh = false;
	bool fillhole = false;
	bool surfaceholefill = false;
	bool avenormal = false;

	bool normalize = false;
	bool remove_hanging_face = false;
	double DeleteDulPoint = -1;
	double reset_face_id_by_angle = -1;
	double scale_factor = 1;
	int reparam_way = -1;
	int shuffle_num = 1;
	int fairing_k = -1;

	vector<int> saveid; /// these ids will not be removed after reset_orient_faceid
	vector<double> rotateVec;
	vector<double> boxVec;
	vector<double> createbox;

	app.add_option("--boundingbox,-b", boxVec, "add a bounding box to a existed mesh. Format is (length ratio, width ratio, hight ratio, refine number), (5,5,5,4) is recommand");
	app.add_option("-r", rotateVec, "input rotate param. Format is (start_x, start_y, start_z, end_x, end_y, end_z, angle) or (end_x, end_y, end_z, angle). angle value scale is (0, 2).");
	app.add_option("-i", input_filenames, "input filename. (string, required, supported format: vtk, mesh, pls, obj)")->required()->expected(1, 3);
	app.add_option("-p", input_filename_ex, "input filename. (string, required)");
	app.add_flag("-k", exportVTK, "Write mesh in VTK format.");
	app.add_flag("-e", exportEpsVTK, "Set eps in VTK format.");
	app.add_flag("-m", exportMESH, "Write mesh in MESH/MEDIT format.");
	app.add_flag("-y", exportPLY, "Write mesh in PLY format.");
	app.add_flag("-s", exportPLS, "Write mesh in PLS format.");
	app.add_flag("-f", exportFacet, "Write mesh in facet format.");
	app.add_flag("-o", exportOBJ, "Write mesh in OBJ format.");
	app.add_flag("--stl", exportSTL, "Write mesh in STL format.");
	app.add_flag("--stl_in", exportStlIn, "Write mesh in stl.in format.");
	app.add_flag("--remesh", bremesh, "Remesh");
	app.add_flag("--fillhole", fillhole, "Fill hole by topology");
	app.add_flag("--surfaceholefill", surfaceholefill, "Fill surfacehole by fairing.");
	app.add_flag("--shuffle", shuffleMark, "Shuffle surface_id for view clearly.");
	app.add_option("--fairing_k", fairing_k, "fairing_k is 2 at least and mustn't be too big.");
	app.add_option("--shuffle_num", shuffle_num, "Shuffle number is [1, 100].");
	app.add_option("--reparam", reparam_way, "input reparameter way. 0 is Tuttle. 1 is harmonic.");
	app.add_option("--create_box", createbox, "Create a mesh only contain a box. Format is (x1_min, y1_min, z1_min, x1_max, y1_max, z1_max, ...)");
	app.add_flag("--reverse_orient", reverseFacetOrient, "Reverse Facet Orient.");
	app.add_flag("--reset_orient", resetOritation, "Regularize orientation");
	app.add_flag("--check_orient", checkOritation, "Check orientation error");
	app.add_flag("--reset_orient_faceid", resetOritationFaceid, "Regularize orientation and reset the facet mask by connected graph component index.");
	app.add_flag("--reset_faceid_angle", reset_face_id_by_angle, "Regularize orientation and reset the facet mask by facet angle.");
	app.add_option("--save_id", saveid, "saved id.");
	app.add_flag("--rm_zero_area", RepairZeroAera, "Repair mesh file for the facet's area that equal to zero.");
	app.add_flag("--rm_box", removebox, "remove the component with the biggest volume. warning, the facet id may be reoriented");
	app.add_flag("--rm_duplicate_point", DeleteDulPoint, "Delete duplicate points; Format is : e.g. --rm_duplicate_point=1e-6.");
	app.add_flag("--normalize", normalize, "normalize the mesh data to [0,1]");
	app.add_flag("--rm_hanging", remove_hanging_face, "remove the hanging face.");
	app.add_flag("--avenormal", avenormal, "cal ave normal.");
	app.add_option("--scale_factor", scale_factor, "Controls the scaling of grid coordinates. A value greater than 1 will enlarge the grid, while a value less than 1 will shrink it.");
	if (resetOritationFaceid)
		resetOritation = false;

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError& e)
	{
		return app.exit(e);
	}

	string input_filename = input_filenames[0];
	size_t input_dotpos = input_filename.find_last_of('.');
	string input_postfix = input_filename.substr(input_dotpos + 1, input_filename.length() - input_dotpos - 1);

	Mesh mesh;
	int cou;
	std::map<int, double> mpd;
	std::map<int, vector<int>> mpi;
	std::string suffix = "";
	Eigen::MatrixXd V_uv;

	if (input_postfix == "vtk")
		MESHIO::readVTK(input_filename, mesh, "surface_id");
	else if (input_postfix == "mesh")
		MESHIO::readMESH(input_filename, mesh);
	else if (input_postfix == "pls")
		MESHIO::readPLS(input_filename, mesh);
	else if (input_postfix == "ply")
		MESHIO::readPLY(input_filename, mesh);
	else if (input_postfix == "obj")
		MESHIO::readOBJ(input_filename, mesh);
	else if (input_postfix == "facet")
		MESHIO::readFacet(input_filename, mesh);
	else if (input_postfix == "stl")
		MESHIO::readSTL(input_filename, mesh);
	else if (input_postfix == "node" || input_postfix == "ele" || input_postfix == "face")
	{
		string nodefilename = "";
		string facefilename = "";
		string elemfilename = "";
		for (auto t_filename : input_filenames)
		{
			size_t t_dotpos = t_filename.find_last_of('.');
			string t_postfix = t_filename.substr(t_dotpos + 1, t_filename.length() - t_dotpos - 1);
			if (t_postfix == "node")
				nodefilename = t_filename;
			else if (t_postfix == "ele")
				elemfilename = t_filename;
			else if (t_postfix == "face")
				facefilename = t_filename;
		}
		if (nodefilename.empty())
		{
			std::cout << "Node file is required." << endl;
		}
		if (!elemfilename.empty())
		{
			MESHIO::readTetgen(nodefilename, elemfilename, mesh);
		}
		else if (!facefilename.empty())
		{
			MESHIO::readTetgen(nodefilename, facefilename, mesh);
		}
	}
	else
	{
		cout << "Unsupported input format - " << input_postfix << endl;
		return -1;
	}

	// cout for checking reading.
	cout << "Read " << mesh.Vertex.rows() << " points." << endl;
	cout << "Read " << mesh.Topo.rows() << " elements." << endl;

	if (exportEpsVTK)
		MESHIO::readEPS(input_filename, cou, mpd, mpi);

	//********* Reparameter *********
	if (reparam_way != -1)
	{
		switch (reparam_way)
		{
		case 0:
			suffix = ".tutte";
			MESHIO::buildTuttleParameter(mesh.Vertex, mesh.Topo, V_uv);
			for (int i = 0; i < V_uv.rows(); i++)
			{
				mesh.Vertex(i, 0) = V_uv(i, 0);
				mesh.Vertex(i, 1) = V_uv(i, 1);
				mesh.Vertex(i, 2) = 0;
			}
			break;
		case 1:
			suffix = ".harmonic";
			MESHIO::buildHarmonicParameter(mesh.Vertex, mesh.Topo, V_uv);
			for (int i = 0; i < V_uv.rows(); i++)
			{
				mesh.Vertex(i, 0) = V_uv(i, 0);
				mesh.Vertex(i, 1) = V_uv(i, 1);
				mesh.Vertex(i, 2) = 0;
			}
			break;
		default:
			break;
		}
	}


	//********* Normalize *********

	if (normalize) {
		MESHIO::Normalize(mesh);
	}


	//********* Ave normal *********
	if (avenormal) {
		auto avenormal = MESHIO::ave_normal(mesh);
		std::cout << "Ave area normal is : " << avenormal.x() << " " << avenormal.y() << " " << avenormal.z() << "\n";
	}


	//********* Rotate *********
	if (!rotateVec.empty())
		MESHIO::rotatePoint(rotateVec, mesh);

	//********* Shuffle **********
	if (shuffleMark)
	{
		std::cout << "shuffle surface start.\n";
		MESHIO::shuffleSurfaceid(shuffle_num, mesh);
		suffix = ".shuffle";
	}

	//********* Add Box *********
	if (!boxVec.empty())
		MESHIO::addBox(boxVec, mesh);

	//********* Regularize mesh oritation *********
	if (resetOritation) {
		auto T = mesh.Topo;
		Eigen::MatrixXi C;
		igl::bfs_orient(T, mesh.Topo, C);
		// MESHIO::resetOrientation(mesh);
	}
	if (resetOritationFaceid)
		MESHIO::resetOrientation(mesh, true, saveid);
	if (checkOritation)
		MESHIO::checkOrientation(mesh);

	//********* modify facet orient ******
	if (reverseFacetOrient)
	{
		MESHIO::reverseOrient(mesh.Topo);
	}

	if (reset_face_id_by_angle > 0) {
		MESHALG::facet_classification(mesh.Vertex, mesh.Topo, reset_face_id_by_angle, mesh.Masks);
	}

	//********* Delete dulplicate point ********
	if (DeleteDulPoint != -1)
	{
		MESHIO::removeDulplicatePoint(mesh.Vertex, mesh.Topo, DeleteDulPoint);
	}

	//********* repair ********
	if (RepairZeroAera)
	{
		MESHIO::repair(mesh);
	}


	if (remove_hanging_face) {
		MESHIO::removeHangingFace(mesh);
	}

	//********* HoleFill *********
	if (fillhole) {
		// MESHIO::topoFillHole(mesh);
		MESHIO::dynamicFillHole(mesh);
	}

	//********* HoleFill *********
	if (surfaceholefill) {
		if (fairing_k < 2)
		{
			std::cout << "Fairing_k is error!Fairing_k at least be 2\n";
			return 0;
		}
		else
		{
			SurfaceHoleFilling holefill(mesh);
			holefill.init_hole_fill();
			//if (fairing_k != 2 && fairing_k != 3)std::cout << "fairing_k must be 2 or 3.\n";
			holefill.fair(fairing_k); // ���ڱ�ַ�
			mesh = holefill.get_all_mesh();
		}


	}

	//********* Remesh *********
	if (bremesh)
	{
		MESHIO::remesh(mesh);
	}



	//********* Create some box ********
	if (createbox.size() != 0)
	{
		if (createbox.size() % 6 != 0)
		{
			std::cout << "Format is error !" << std::endl;
			return -1;
		}
		MESHIO::createBox(createbox, mesh);
	}

	//********* remove box ********
	if (removebox) {
		MESHIO::removeBox(mesh);
	}
	//********* export ********

	if (std::abs(scale_factor - 1) > 1e-8) {
		MESHIO::scale(mesh, scale_factor);
	}

	if (exportVTK)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.vtk";
		MESHIO::writeVTK(output_filename, mesh, "surface_id");
	}
	if (exportMESH)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.mesh";
		MESHIO::writeMESH(output_filename, mesh);
	}
	if (exportPLY)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.ply";
		MESHIO::writePLY(output_filename, mesh);
	}
	if (exportPLS)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.pls";
		MESHIO::writePLS(output_filename, mesh);
	}
	if (exportFacet)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.facet";
		MESHIO::writeFacet(output_filename, mesh);
	}
	if (exportEpsVTK)
	{ // This is to generate AutoGrid to control local eps.
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".eps.vtk";
		MESHIO::writeEpsVTK(output_filename, mesh, cou, mpd, mpi);
	}
	if (exportOBJ)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.obj";
		MESHIO::writeOBJ(output_filename, mesh);
	}
	if (exportStlIn)
	{
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + "_stl.in";
		MESHIO::writeStlIn(output_filename, mesh);
	}
	if (exportSTL)
	{	
		string output_filename = input_filename.substr(0, input_dotpos) + suffix + ".o.stl";
		MESHIO::writeSTL(output_filename, mesh);
	}

	cout << "Write " << mesh.Vertex.rows() << " points." << endl;
	cout << "Write " << mesh.Topo.rows() << " elements." << endl;

	return 0;
}

