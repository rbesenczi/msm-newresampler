/*
Copyright (c) 2022 King's College London, MeTrICS Lab, Renato Besenczi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include "mesh.h"

using namespace std;

namespace newresampler {

Mesh::Mesh() {
    string default_space = "NIFTI_XFORM_TALAIRACH";
    std::vector<double> transform(16, 0);
    transform[0] = 1;
    transform[5] = 1;
    transform[10] = 1;
    transform[15] = 1;
    NEWMESH::GIFTIcoordinateSystem default_coord(default_space, default_space, transform); // save out a default coord system
    global_defaultcoord.push_back(default_coord);
}

Mesh::Mesh(const Mesh &m) : normals(m.normals), pvalues(m.pvalues), tvalues(m.tvalues),
                            global_metaData(m.global_metaData), global_Attributes(m.global_Attributes),
                            global_GIFTIlabels(m.global_GIFTIlabels), global_defaultcoord(m.global_defaultcoord){

    for (const auto & point : m.points) {
        std::shared_ptr<Mpoint> pt = std::make_shared<Mpoint>(*point);
        points.push_back(pt);
    }

    for (const auto& triangle : m.triangles) {
        int v0 = triangle.get_vertex_no(0),
            v1 = triangle.get_vertex_no(1),
            v2 = triangle.get_vertex_no(2);
        Triangle tr(points[v0], points[v1], points[v2], triangle.get_no());
        triangles.push_back(tr);
    }

    for (int i = 0; i < (int)points.size(); i++)
        normals.push_back(local_normal(i));
}

Mesh &Mesh::operator=(const Mesh &m) {

    points.clear();
    triangles.clear();
    normals = m.normals;
    pvalues = m.pvalues;
    tvalues = m.tvalues;

    global_metaData = m.global_metaData;
    global_Attributes = m.global_Attributes;
    global_GIFTIlabels = m.global_GIFTIlabels;
    global_defaultcoord = m.global_defaultcoord;

    for (const auto& point : m.points) {
        std::shared_ptr<Mpoint> pt = std::make_shared<Mpoint>(*point);
        points.push_back(pt);
    }

    for (const auto & triangle : m.triangles) {
        int v0 = triangle.get_vertex_no(0), v1 = triangle.get_vertex_no(1), v2 = triangle.get_vertex_no(2);
        Triangle tr(points[v0], points[v1], points[v2], triangle.get_no());
        triangles.push_back(tr);
    }

    return *this;
}

Mesh::Mesh(Mesh&& m) noexcept {

    points = std::move(m.points);
    triangles = std::move(m.triangles);
    normals = std::move(m.normals);
    pvalues = std::move(m.pvalues);
    tvalues = std::move(m.tvalues);

    global_metaData = std::move(m.global_metaData);
    global_Attributes = std::move(m.global_Attributes);
    global_GIFTIlabels = std::move(m.global_GIFTIlabels);
    global_defaultcoord =  std::move(m.global_defaultcoord);
}

Mesh& Mesh::operator=(Mesh&& m) noexcept {

    points.clear();
    triangles.clear();

    points = std::move(m.points);
    triangles = std::move(m.triangles);
    normals = std::move(m.normals);
    pvalues = std::move(m.pvalues);
    tvalues = std::move(m.tvalues);

    global_metaData = std::move(m.global_metaData);
    global_Attributes = std::move(m.global_Attributes);
    global_GIFTIlabels = std::move(m.global_GIFTIlabels);
    global_defaultcoord =  std::move(m.global_defaultcoord);

    return *this;
}

void Mesh::push_triangle(const Triangle &t) {

    std::vector<int> _nos;
    triangles.push_back(t);

    for (int i = 0; i < 3; i++) {
        points[t.get_vertex(i).get_no()]->push_triangle(t.get_no());
        _nos.push_back(t.get_vertex(i).get_no());
    }

    if (!points[_nos[0]]->is_neighbour(_nos[1]))
        points[_nos[0]]->push_neighbour(_nos[1]); // added by emma to ensure vertex neighbourhoods are saved
    if (!points[_nos[0]]->is_neighbour(_nos[2])) points[_nos[0]]->push_neighbour(_nos[2]);

    if (!points[_nos[1]]->is_neighbour(_nos[0])) points[_nos[1]]->push_neighbour(_nos[0]);
    if (!points[_nos[1]]->is_neighbour(_nos[2])) points[_nos[1]]->push_neighbour(_nos[2]);

    if (!points[_nos[2]]->is_neighbour(_nos[0])) points[_nos[2]]->push_neighbour(_nos[0]);
    if (!points[_nos[2]]->is_neighbour(_nos[1])) points[_nos[2]]->push_neighbour(_nos[1]);
}

Point Mesh::local_normal(int pt) const {

    Point v;
    for (int i = 0; i < points[pt]->ntriangles(); i++)
        v += triangles[points[pt]->get_trID(i)].normal();
    v.normalize();

    return v;
}

void Mesh::estimate_normals() {

    normals.clear();
    normals.reserve(points.size());

    for (int i = 0; i < (int)points.size(); i++)
        normals.emplace_back(local_normal(i));
}

vector<float> Mesh::getPointsAsVectors() const {

    vector<float> ret;
    for (int i = 0; i < nvertices(); i++)
    {
        ret.emplace_back(get_point(i).get_coord().X);
        ret.emplace_back(get_point(i).get_coord().Y);
        ret.emplace_back(get_point(i).get_coord().Z);
    }
    return ret;
}

vector<float> Mesh::getValuesAsVectors(int d) const {

    vector<float> ret;
    for (unsigned int i = 0; i < pvalues[d].size(); i++)
        ret.push_back(get_pvalue(i, d));

    return ret;
}

vector<int> Mesh::getTrianglesAsVector() const {

    vector<int> ret;
    for (int i = 0; i < ntriangles(); i++)
    {
        ret.emplace_back(triangles[i].get_vertex(0).get_no());
        ret.emplace_back(triangles[i].get_vertex(1).get_no());
        ret.emplace_back(triangles[i].get_vertex(2).get_no());
    }
    return ret;
}

void Mesh::set_pvalue(unsigned int i, float val, int dim) {

    if (pvalues.empty()) {
        cout << "Mesh::set_pvalue, warning, pvalues have not been initialised. Generating pvalues field of length=number of vertices" << endl;
        initialize_pvalues(dim + 1, true);
    } else if ((int) pvalues.size() < dim + 1) {
        cout << "Mesh::set_pvalue, warning, index exceeds known pvalues dimension. Appending pvalues" << endl;
        initialize_pvalues(dim - pvalues.size() + 1, true);
    }
    if (i >= pvalues[dim].size()) {
        cout << i << " dim " << dim << " " << pvalues[dim].size() << endl;
        throw MeshException("Mesh::set_pvalue, index is incompatible with data dimensions");
    }

        pvalues[dim][i] = val;
} // this will be problematic is i is > the size of the vector

void Mesh::set_pvalues(const NEWMAT::Matrix &M, bool appendFieldData) {
    bool verticesAreColumns = false;

    if (!points.empty()) { // check that if points have been supplied, that the data has the same length
        if (M.Ncols() == (int) points.size()) verticesAreColumns = true;
        else if (M.Nrows() == (int) points.size()) verticesAreColumns = false;
        else throw MeshException(" Cannot assign data to mesh. Dimensions do not match that of mesh");
    } else if (!pvalues.empty()) {
        if (M.Ncols() == (int) pvalues[0].size()) verticesAreColumns = true;
        else if (M.Nrows() == (int) pvalues[0].size()) verticesAreColumns = false;
        else throw MeshException(" Cannot assign data to mesh. Dimensions do not match that of mesh");
    } else if (M.Ncols() > M.Nrows()) verticesAreColumns = true;

    if (!appendFieldData)
        pvalues.clear();

    if (verticesAreColumns) {
        for (int i = 1; i <= M.Nrows(); i++) {
            vector<float> tmp_pvalues;
            for (int j = 1; j <= M.Ncols(); j++)
                tmp_pvalues.push_back(M(i, j));

            pvalues.push_back(tmp_pvalues);
        }
    } else {
        for (int i = 1; i <= M.Ncols(); i++) {
            vector<float> tmp_pvalues;
            for (int j = 1; j <= M.Nrows(); j++)
                tmp_pvalues.push_back(M(j, i));

            pvalues.push_back(tmp_pvalues);
        }
    }
}

void Mesh::initialize_pvalues(int dim, bool appendFieldData) {

    if (!appendFieldData)
        pvalues.clear();

    pvalues.resize(dim + pvalues.size(), std::vector<float>(points.size(), 0.0));
}

double Mesh::calculate_MaxVD() const {

    double dist, MaxVD = std::numeric_limits<double>::lowest();
    const double RAD = 100.0;

    for (long unsigned int i = 0; i < points.size(); i++)
    {
        Point control_point = points[i]->get_coord();
        for (auto it = nbegin(i); it != nend(i); it++)
        {
            dist = 2 * RAD * asin((control_point - get_coord(*it)).norm() / (2 * RAD));
            if(dist > MaxVD) MaxVD = dist;
        }
    }
    return MaxVD;
}

double Mesh::Calculate_MVD() const {

    // mean_inter_vertex_distance
    // calculate from first vertex only so not true unless regularised!
    double k = 0.0, kr = 0.0;
    for (int i = 0; i < nvertices(); i++)
    {
        Point Sp = get_coord(i);
        for (auto j = nbegin(i); j != nend(i); j++)
        {
            k++;
            kr += (get_coord(*j) - Sp).norm();
        }
    }

    return kr / k;
}

Mesh::FileType Mesh::meshFileType(const string &filename) const {
    // 1: ascii, 2:vtk, 3: gii, 4 as .txt (for supplyng data in a text file) -1: unknown
    if (filename.size() <= 5) { return FileType::DEFAULT; }
    string last_3 = filename.substr(filename.size() - 3, 3);
    if (last_3 == ".gz") {
        last_3 = filename.substr(filename.size() - 6, 3);
        if (last_3 == "gii") { return FileType::GIFTI; }
    }

    if (last_3 == "gii") { return FileType::GIFTI; }
    else if (last_3 == "txt") { return FileType::MATRIX; }
    else if (last_3 == "dpv") { return FileType::DPV; }
    else if (last_3 == "asc") { return FileType::ASCII; }
    else if (last_3 == "vtk") { return FileType::VTK; }

    ifstream f(filename.c_str());
    //reading the header
    string header;
    getline(f, header);
    {
        string::size_type pos = header.find("# vtk DataFile Version");
        if (pos != string::npos) {
            f.close();
            return FileType::VTK;
        }
    }
    {
        string::size_type pos = header.find("#!ascii");
        if (pos != string::npos) {
            f.close();
            return FileType::ASCII;
        }
    }
    return FileType::DEFAULT;
}

void Mesh::load(const string &filename, const bool loadSurfaceData, const bool appendFieldData) {
    //  appendFieldData refers only to pvalue/tvalue fields.
    //  As such if you are reading a .func or a .shape following a .surf appendFieldData=false.
    Mesh::FileType type = meshFileType(filename);
    if (type == FileType::ASCII) {
        load_ascii(filename, loadSurfaceData, appendFieldData);
    } else if (type == FileType::VTK) {
        load_vtk(filename);
    } else if (type == FileType::GIFTI) {
        load_gifti(filename, loadSurfaceData, appendFieldData);
    } else if (type == FileType::MATRIX || type == FileType::DPV) {
        load_matrix(filename, type);
    } else {
        cerr << "Mesh::load:error reading file: " << filename << "  ... Unknown format" << endl;
        exit(1);
    }
}

void Mesh::load_gifti(const string &filename, const bool loadSurfaceData, const bool appendFieldData) {
    NEWMESH::GIFTIwrapper reader;
    reader.readGIFTI(filename);

    if (!appendFieldData) {
        pvalues.clear();
        tvalues.clear();
    }

    if (loadSurfaceData) {
        points.clear();
        triangles.clear();
        vector<NEWMESH::GIFTIfield> surfaceData = reader.returnSurfaceFields();
        global_defaultcoord = surfaceData[0].getCoordSystems();  // save out a default coord system

        for (int point = 0; point < surfaceData[0].getDim(0); point++) {
            vector<float> coords(surfaceData[0].fVector(point));
            std::shared_ptr<Mpoint> m = std::make_shared<Mpoint>(
                    coords[0], coords[1], coords[2], point);
            points.push_back(m);
        }

        for (int triangle = 0; triangle < surfaceData[1].getDim(0); triangle++) {
            vector<int> points_vec(surfaceData[1].iVector(triangle));
            Triangle t(points[points_vec[0]], points[points_vec[1]], points[points_vec[2]], triangle);
            push_triangle(t);
        }
    }

    vector<NEWMESH::GIFTIfield> nonSurfaceData = reader.returnNonSurfaceFields();
    int point;

    for (auto & dim : nonSurfaceData) {
        if (points.empty() || dim.getDim(0) == (int) points.size()) {
            vector<float> tmp_pvalues;
            for (point = 0; point < dim.getDim(0); point++) {
                tmp_pvalues.push_back(dim.fScalar(point));
            }
            pvalues.push_back(tmp_pvalues);
        } else if (dim.getDim(0) ==
                   (int) triangles.size()) {
            vector<float> tmp_tvalues;
            for (point = 0; point < dim.getDim(0); point++) {
                tmp_tvalues.push_back(dim.fScalar(point));
            }
        } else { throw MeshException(" mismatch between data and surface dimensions"); }
    }
    global_metaData = reader.metaData;
    global_Attributes = reader.extraAttributes; //Any other attributes that the tag has
    global_GIFTIlabels = reader.GIFTIlabels;
}

void Mesh::load_vtk(const string &filename) {
    // cannot provide field data - surface only
    points.clear();
    triangles.clear();

    ifstream f(filename.c_str());
    if (f.is_open()) {
        //reading the header
        string header;
        getline(f, header);
        string::size_type pos = header.find("# vtk DataFile Version");
        if (pos == string::npos) {
            cerr << "Mesh::load_vtk:error in the header" << endl;
            exit(1);
        }
        getline(f, header);
        getline(f, header);
        getline(f, header);
        int NVertices, NFaces;
        f >> header >> NVertices >> header;
        //reading the points
        for (int i = 0; i < NVertices; i++) {
            double x, y, z;
            f >> x >> y >> z;
            std::shared_ptr<Mpoint> m = std::make_shared<Mpoint>(x, y, z, i);
            points.push_back(m);
        }
        f >> header >> NFaces >> header;

        //reading the triangles
        for (int i = 0; i < NFaces; i++) {
            int p0, p1, p2;
            int j;
            f >> j >> p0 >> p1 >> p2;
            Triangle t(points[p0], points[p1], points[p2], i);
            push_triangle(t);

        }
        f >> header >> header;
        f >> header >> header >> header;
        f >> header >> header;
        //reading the values

        for (int i = 0; i < NVertices; i++) {
            int val;
            f >> val;
            set_pvalue(i, val);
        }
        f.close();
    } else {
        cout << "Mesh::error opening file: " << filename << endl;
        exit(1);
    }
}

void Mesh::load_ascii(const string &filename, const bool loadSurfaceData,
                      const bool appendFieldData) { //load a freesurfer ascii mesh
// if loadSurfaceData = true then pvalues are ignored (for compatability with gifti, where .surf and .func are loaded separately)

    if (loadSurfaceData) {
        triangles.clear();
        points.clear();
    }

    if (!appendFieldData) {
        pvalues.clear();
        tvalues.clear();
    }

    ifstream f(filename.c_str());
    if (f.is_open()) {
        //reading the header
        string header;
        getline(f, header);
        string::size_type pos = header.find("#!ascii");
        if (pos == string::npos) {
            cerr << "Mesh::load_ascii:error in the header" << endl;
            exit(1);
        }

        //reading the size of the mesh
        int NVertices, NFaces;
        f >> NVertices >> NFaces;

        if (!loadSurfaceData) {
            vector<float> tmp_pvalues(NVertices, 0);
            pvalues.push_back(tmp_pvalues);
        }
        //reading the points
        for (int i = 0; i < NVertices; i++) {
            double x, y, z;
            float val;
            f >> x >> y >> z >> val;

            if (loadSurfaceData) {
                std::shared_ptr<Mpoint> m = std::make_shared<Mpoint>(x, y, z, i);
                points.push_back(m);
            } else pvalues[pvalues.size() - 1][i] = val;
        }
        //reading the triangles

        if (!loadSurfaceData) {
            vector<float> tmp_tvalues(NFaces, 0);
            tvalues.push_back(tmp_tvalues);
        }

        for (int i = 0; i < NFaces; i++) {
            int p0, p1, p2;
            float val;
            f >> p0 >> p1 >> p2 >> val;

            if (loadSurfaceData) {
                Triangle t(points[p0], points[p1], points[p2], i);
                push_triangle(t);
            } else tvalues[tvalues.size() - 1][i] = val;
        }
        f.close();
    } else {
        cout << "Mesh::load_ascii:error opening file: " << filename << endl;
        exit(1);
    }
}

void Mesh::load_ascii_file(const string &filename) { //load a freesurfer ascii mesh for pvalues only
    // loads surface & data
    clear();

    ifstream f(filename.c_str());
    if (f.is_open()) {
        //reading the header
        string header;
        getline(f, header);
        string::size_type pos = header.find("#!ascii");
        if (pos == string::npos) {
            cerr << "Mesh::load_ascii:error in the header" << endl;
            exit(1);
        }

        //reading the size of the mesh
        int NVertices, NFaces;
        f >> NVertices >> NFaces;

        vector<float> tmp_pvalues;

        for (int i = 0; i < NVertices; i++) {
            double x, y, z;
            float val;
            f >> x >> y >> z >> val;
            std::shared_ptr<Mpoint> m = std::make_shared<Mpoint>(x, y, z, i);
            points.push_back(m);
            tmp_pvalues.push_back(val);
        }
        pvalues.push_back(tmp_pvalues);

        vector<float> tmp_tvalues;
        tvalues.push_back(tmp_tvalues);

        for (int i = 0; i < NFaces; i++) {
            int p0, p1, p2;
            float val;
            f >> p0 >> p1 >> p2 >> val;
            Triangle t(points[p0], points[p1], points[p2], i);
            push_triangle(t);
            tmp_tvalues.push_back(val);
        }
        tvalues.push_back(tmp_tvalues);
        f.close();

    } else {
        cout << "Mesh::load_ascii:error opening file: " << filename << endl;
        exit(1);
    }
}

void Mesh::load_matrix(const string &filename,
                       const Mesh::FileType &type) {  // for pvalues only - when data is held in a textfile
    // cannot provide field data - surface only (also reads .dpv files)
    NEWMAT::Matrix tmp, tmp2;
    tmp = MISCMATHS::read_ascii_matrix(filename);
    if (type == FileType::DPV) {
        tmp2 = tmp;
        tmp.ReSize(tmp.Nrows(), 1);

        if (tmp2.Ncols() != 5) {
            cout << "Mesh::load_dpv:error opening file (wrong format) : " << filename << endl;
            exit(1);
        }
        for (int i = 1; i <= tmp.Nrows(); i++) {
            if (tmp2(i, 1) != i - 1) {
                cout << i << " " << tmp2(i, 1) << "Mesh::load_dpv:error opening file (wrong format 2) : "
                     << filename << endl;
                exit(1);
            }
            tmp(i, 1) = tmp2(i, 5);
        }
    }
    if (tmp.Nrows() == 0 && tmp.Ncols() == 0) {
        cout
                << "Mesh::load_txt:error opening file (wrong format). Note matrix file must be delimited with spaces ?? "
                << endl;
        exit(1);
    }
    if (tmp.Ncols() == 5) {
        int ind = 0;
        for (int i = 1; i <= tmp.Nrows(); i++)
            if (tmp(i, 1) == i - 1) ind++;

        if (ind == tmp.Nrows()) {
            cout << "WARNING: this looks like a dpv file but is being read as a text file!! " << endl;
        }
    }
    set_pvalues(tmp);
}

void Mesh::save(const string &filename) const {

    Mesh::FileType type = meshFileType(filename);
    switch (type) {
        case Mesh::FileType::DPV:
            save_dpv(filename);
            break;
        case Mesh::FileType::MATRIX:
            save_matrix(filename);
            break;
        case Mesh::FileType::ASCII:
            save_ascii(filename);
            break;
        case Mesh::FileType::VTK:
            save_vtk(filename);
            break;
        case Mesh::FileType::DEFAULT:
            save_gifti(filename);
            break;
        default:
            throw MeshException {"invalid file type "};
    }
}

void Mesh::save_gifti(const string &s) const {
    string filename(s);
    string last_3 = filename.substr(filename.size() - 3, 3);
    NEWMESH::GIFTIwrapper writer;

    string subtype;

    if (last_3 != "gii") {
        if (last_3 != ".gz") {
            filename = filename + ".gii";
        }
        subtype = filename.substr(filename.size() - 9, 5);

    } else subtype = filename.substr(filename.size() - 9, 5);

    writer.metaData = global_metaData;
    writer.extraAttributes = global_Attributes; //Any other attributes that the tag has
    writer.GIFTIlabels = global_GIFTIlabels;
    //writer.report();
    if (subtype == ".surf") {

        vector<NEWMESH::GIFTIfield> surfaceData;
        vector<int> dims(2, 3); //default as 3-vector for surface
        dims[0] = points.size();

        vector<float> pointData(getPointsAsVectors());
        surfaceData.emplace_back(NiftiIO::NIFTI_INTENT_POINTSET, NiftiIO::NIFTI_TYPE_FLOAT32, 2, dims.data(), pointData.data(),
                           GIFTI_IND_ORD_ROW_MAJOR, global_defaultcoord);

        dims[0] = triangles.size();
        vector<int> triangleData(getTrianglesAsVector());

        surfaceData.emplace_back(NiftiIO::NIFTI_INTENT_TRIANGLE, NiftiIO::NIFTI_TYPE_INT32, 2, dims.data(), triangleData.data(),
                           GIFTI_IND_ORD_ROW_MAJOR, global_defaultcoord);
        writer.allFields = surfaceData;
        writer.writeGIFTI(filename, GIFTI_ENCODING_B64GZ);
    } else if (subtype == ".func" || subtype == "shape") {
        vector<int> scalardims(1);
        vector<NEWMESH::GIFTIfield> fieldData;

        for (unsigned int dim = 0; dim < pvalues.size(); dim++) {
            scalardims[0] = pvalues[dim].size();
            vector<float> valueData(getValuesAsVectors(dim));
            fieldData.emplace_back(NiftiIO::NIFTI_INTENT_NONE, NiftiIO::NIFTI_TYPE_FLOAT32, 1, scalardims.data(), valueData.data(),
                               GIFTI_IND_ORD_ROW_MAJOR);
        }
        writer.allFields = fieldData;
        writer.writeGIFTI(filename, GIFTI_ENCODING_B64GZ);
    }
}

void Mesh::save_vtk(const string &s) const {
    string filename(s);
    string last_3 = filename.substr(filename.size() - 3, 3);
    if (last_3 != "vtk") {
        filename = filename + ".vtk";
    }

    ofstream flot(filename.c_str());
    if (flot.is_open()) {
        flot << "# vtk DataFile Version 3.0" << endl
             << "surface file" << endl
             << "ASCII" << endl
             << "DATASET POLYDATA" << endl
             << "POINTS ";
        flot << points.size() << "  float" << endl;

        for (const auto & point : points) {
            //	flot.precision(6);
            flot << point->get_coord().X << " "
                 << point->get_coord().Y << " "
                 << point->get_coord().Z << endl;
#ifdef PPC64
            if ((n++ % 20) == 0) flot.flush();
#endif
        }
        flot << "POLYGONS " << triangles.size() << " " << triangles.size() * 4 << endl;
        for (const auto & triangle : triangles)
            flot << "3 "
                 << triangle.get_vertex(0).get_no() << " "
                 << triangle.get_vertex(1).get_no() << " "
                 << triangle.get_vertex(2).get_no() << " " << endl;
#ifdef PPC64
        if ((n++ % 20) == 0) flot.flush();
#endif
    } else {
        cerr << "::save_vtk:error opening file " << filename << " for writing" << endl;
        exit(1);
    }
}

void Mesh::save_ascii(const string &s) const {
    string filename(s);
    string last_3 = filename.substr(filename.size() - 3, 3);
    if (last_3 != "asc") { filename = filename + ".asc"; }

    ofstream f(filename.c_str());
    stringstream flot;
    if (f.is_open()) {
        int ptcount(0), tricount(0);
        for (unsigned int i = 0; i < points.size(); i++) {
            float val;
            if (pvalues.empty()) val = 0;
            else val = pvalues[0][i]; // value of first column

            flot << points[i]->get_coord().X << " "
                 << points[i]->get_coord().Y << " "
                 << points[i]->get_coord().Z << " "
                 << val << endl;
            ptcount++;
        }
        for (unsigned int i = 0; i < triangles.size(); i++) {
            float val;
            if (tvalues.empty()) val = 0;
            else val = tvalues[0][i]; // value of first column

            flot << triangles[i].get_vertex(0).get_no() << " "
                 << triangles[i].get_vertex(1).get_no() << " "
                 << triangles[i].get_vertex(2).get_no() << " " << val << endl;
            tricount++;
        }
        f << "#!ascii from Mesh" << endl;
        f << ptcount << " " << tricount << endl << flot.str();
        f.close();
    } else cerr << "Mesh::save_ascii:error opening file for writing: " << s << endl;

}

void Mesh::save_dpv(const string &s) const {
    string filename(s);
    string last_3 = filename.substr(filename.size() - 3, 3);
    if (last_3 != "dpv") { filename = filename + ".dpv"; }

    ofstream f(filename.c_str());
    if (f.is_open()) {
        if (pvalues.empty()) {
            throw MeshException("Mesh::save_dpv, cannot write out as dpv as there is no data");
        } else {
            if (points.size() != pvalues[0].size()) {
                throw MeshException("Mesh::save_dpv, data and mesh dimensions do not agree");
            }

            for (unsigned int i = 0; i < pvalues[0].size(); i++) {
                float val = pvalues[0][i];  /// outputs value for first column only
                if (i < 100)
                    f << setfill('0') << setw(3) << i;
                else
                    f << i;

                if (!points.empty()) {
                    f << " " << points[i]->get_coord().X
                      << " " << points[i]->get_coord().Y
                      << " " << points[i]->get_coord().Z
                      << " " << val << endl;

                } else {
                    f << " 0 0 0 " << val << endl;

                }
            }
            f.close();
        }
    } else cerr << "Mesh::save_ascii:error opening file for writing: " << s << endl;
}

void Mesh::save_matrix(const string &s) const {
    string filename(s);
    string last_3 = filename.substr(filename.size() - 3, 3);
    if (last_3 != "txt") { filename = filename + ".txt"; }

    ofstream f(filename.c_str());
    if (f.is_open()) {
        if (pvalues.empty()) {
            throw MeshException("Mesh::save_matrix, cannot write out matrix as there is no data");
        } else {
            if (points.size() != pvalues[0].size()) {
                throw MeshException("Mesh::save_matrix, data and mesh dimensions do not agree");
            }
            for (const auto & pvalue : pvalues) {
                for (float i : pvalue) {
                    f << i << " ";
                }
                f << endl;
            }
            f.close();
        }
    } else cerr << "Mesh::save_ascii:error opening file for writing: " << s << endl;
}

void Mesh::clear() {
    points.clear();
    triangles.clear();
    pvalues.clear();
    tvalues.clear();
}

void Mesh::clear_data() {
    pvalues.clear();
    tvalues.clear();
}

const Point &Mesh::get_coord(int n) const {
    if (n >= (int) points.size() || n < 0 || points.empty())
        throw MeshException("get_coord: index exceeds data dimensions");
    else return points[n]->get_coord();
}

float Mesh::get_pvalue(int i, int dim) const {
    if (dim >= (int) pvalues.size() || (int) pvalues[dim].size() < i)
        throw MeshException("get_pvalue: index exceeds data dimensions");
    return pvalues[dim][i];
}

int Mesh::get_total_triangles(int i) const {
    if (i >= (int) points.size() || points.empty())
        throw MeshException("get_coord: index exceeds data dimensions");
    return points[i]->ntriangles();
}

std::vector<Point> Mesh::get_bounding_box() const {

    const double LOWEST = std::numeric_limits<double>::lowest();
    const double LARGEST = std::numeric_limits<double>::max();

    std::vector<Point> bounding_box;
    Point min(LARGEST,LARGEST,LARGEST), max(LOWEST,LOWEST,LOWEST);

    for(const auto& tr : triangles)
    {
        for(int i = 0; i < 3; ++i)
        {
            Point vertex = tr.get_vertex_coord(i);
            if(vertex.X < min.X) min.X = vertex.X;
            if(vertex.Y < min.Y) min.Y = vertex.Y;
            if(vertex.Z < min.Z) min.Z = vertex.Z;
            if(vertex.X > max.X) max.X = vertex.X;
            if(vertex.Y > max.Y) max.Y = vertex.Y;
            if(vertex.Z > max.Z) max.Z = vertex.Z;
        }
    }

    bounding_box.push_back(min);
    bounding_box.push_back(max);

    return bounding_box;
}

int Mesh::get_resolution() const {

    switch(points.size())
    {
        case 42:
            return 1;
        case 162:
            return 2;
        case 642:
            return 3;
        case 2562:
            return 4;
        case 10242:
            return 5;
        case 40962:
            return 6;
        default:
            throw MeshException(
                    " Mesh::get_resolution mesh is not an icosphere or has not been initialised");
    }
}

Point Mesh::estimate_origin() const {

    std::vector<Point> p(4);

    double m11, m12, m13, m14;

    for (int i = 1; i <= 4; i++)
        p[i - 1] = points[std::floor(points.size() / i) - 1]->get_coord();

    NEWMAT::Matrix a(4,4);
    Point c;

    //equation of a sphere as determinant of its variables and 4 sampled points
    for (int i = 0; i < 4; i++)
    {   //find minor 11
        a(i+1,1) = p[i].X;
        a(i+1,2) = p[i].Y;
        a(i+1,3) = p[i].Z;
        a(i+1,4) = 1;
    }
    m11 = a.Determinant();

    for (int i = 0; i < 4; i++)
    {   //find minor 12
        a(i+1,1) = p[i].X * p[i].X + p[i].Y * p[i].Y + p[i].Z * p[i].Z;
        a(i+1,2) = p[i].Y;
        a(i+1,3) = p[i].Z;
        a(i+1,4) = 1;
    }
    m12 = a.Determinant();

    for (int i = 0; i < 4; i++)
    {   //find minor 13
        a(i+1,1) = p[i].X * p[i].X + p[i].Y * p[i].Y + p[i].Z * p[i].Z;
        a(i+1,2) = p[i].X;
        a(i+1,3) = p[i].Z;
        a(i+1,4) = 1;
    }
    m13 = a.Determinant();

    for (int i = 0; i < 4; i++)
    {   //find minor 14
        a(i+1,1) = p[i].X * p[i].X + p[i].Y * p[i].Y + p[i].Z * p[i].Z;
        a(i+1,2) = p[i].X;
        a(i+1,3) = p[i].Y;
        a(i+1,4) = 1;
    }
    m14  = a.Determinant();

    for (int i = 0; i < 4; i++)
    {   // find minor 15
        a(i+1,1) = p[i].X * p[i].X + p[i].Y * p[i].Y + p[i].Z * p[i].Z;
        a(i+1,2) = p[i].X;
        a(i+1,3) = p[i].Y;
        a(i+1,4) = p[i].Y;
    }

    if (m11 != 0.0)
    {
        c.X =  0.5 * (m12 / m11); //center of sphere
        c.Y = -0.5 * (m13 / m11);
        c.Z =  0.5 * (m14 / m11);
    }

    return c;
}

void retessellate(Mesh &mesh) {

    std::vector<std::shared_ptr<Mpoint>> added_points;
    std::vector<Triangle> tr = mesh.get_all_triangles();

    int count = 0;
    int tot_triangles = 0;

    mesh.clear_triangles();

    for (auto &_point: mesh.get_all_points())
        _point->clear();

    for (auto &t: tr) {

        std::shared_ptr<Mpoint> v0 = t.get_vertex_ptr(0);
        std::shared_ptr<Mpoint> v1 = t.get_vertex_ptr(1);
        std::shared_ptr<Mpoint> v2 = t.get_vertex_ptr(2);

        // creates new points at the center of existing faces
        Point pt0((v1->get_coord().X + v2->get_coord().X) / 2,
                  (v1->get_coord().Y + v2->get_coord().Y) / 2,
                  (v1->get_coord().Z + v2->get_coord().Z) / 2);
        Point pt2((v0->get_coord().X + v1->get_coord().X) / 2,
                  (v0->get_coord().Y + v1->get_coord().Y) / 2,
                  (v0->get_coord().Z + v1->get_coord().Z) / 2);
        Point pt1((v0->get_coord().X + v2->get_coord().X) / 2,
                  (v0->get_coord().Y + v2->get_coord().Y) / 2,
                  (v0->get_coord().Z + v2->get_coord().Z) / 2);

        std::shared_ptr<Mpoint> p1, p2, p0;

        bool b0 = true, b1 = true, b2 = true;
        count = 0;
        int index = 0;
        for (const auto &added_point: added_points) {
            index++;
            Point current = added_point->get_coord();
            if (pt0 == current) {
                b0 = false;
                p0 = added_point;
            }
            if (pt1 == current) {
                b1 = false;
                p1 = added_point;
            }
            if (pt2 == current) {
                b2 = false;
                p2 = added_point;
            }
        }

        if (b0) {
            p0 = std::make_shared<Mpoint>(pt0, mesh.nvertices() + count);
            count++;
        }
        if (b1) {
            p1 = std::make_shared<Mpoint>(pt1, mesh.nvertices() + count);
            count++;
        }
        if (b2) {
            p2 = std::make_shared<Mpoint>(pt2, mesh.nvertices() + count);
            count++;
        }

        if (b0) {
            mesh.push_point(p0);
            added_points.push_back(p0);
        }
        if (b1) {
            mesh.push_point(p1);
            added_points.push_back(p1);
        }
        if (b2) {
            mesh.push_point((p2));
            added_points.push_back(p2);
        }

        Triangle t0(p2, p0, p1, tot_triangles);
        tot_triangles++;
        Triangle t1(p1, v0, p2, tot_triangles);
        tot_triangles++;
        Triangle t2(p0, v2, p1, tot_triangles);
        tot_triangles++;
        Triangle t3(p2, v1, p0, tot_triangles);
        tot_triangles++;

        mesh.push_triangle(t0);
        mesh.push_triangle(t1);
        mesh.push_triangle(t2);
        mesh.push_triangle(t3);
    }

    for (auto i = mesh.vbegin(); i != mesh.vend(); i++)
        (*i)->normalize();
}

void retessellate(Mesh& mesh, std::vector<std::vector<int>>& old_tr_nbours) {

    std::vector<std::shared_ptr<Mpoint>> added_points;
    std::vector<Triangle> tr = mesh.get_all_triangles();

    old_tr_nbours.clear();
    old_tr_nbours.reserve(tr.size());

    int count = 0;
    int tot_triangles = 0;

    mesh.clear_triangles();

    for (auto &_point: mesh.get_all_points())
        _point->clear();

    for (auto &t: tr) {

        std::shared_ptr<Mpoint> v0 = t.get_vertex_ptr(0);
        std::shared_ptr<Mpoint> v1 = t.get_vertex_ptr(1);
        std::shared_ptr<Mpoint> v2 = t.get_vertex_ptr(2);

        // creates new points at the center of existing faces
        Point pt0((v1->get_coord().X + v2->get_coord().X) / 2,
                  (v1->get_coord().Y + v2->get_coord().Y) / 2,
                  (v1->get_coord().Z + v2->get_coord().Z) / 2);
        Point pt2((v0->get_coord().X + v1->get_coord().X) / 2,
                  (v0->get_coord().Y + v1->get_coord().Y) / 2,
                  (v0->get_coord().Z + v1->get_coord().Z) / 2);
        Point pt1((v0->get_coord().X + v2->get_coord().X) / 2,
                  (v0->get_coord().Y + v2->get_coord().Y) / 2,
                  (v0->get_coord().Z + v2->get_coord().Z) / 2);

        std::shared_ptr<Mpoint> p1, p2, p0;

        bool b0 = true, b1 = true, b2 = true;
        count = 0;
        int index = 0;
        for (const auto &added_point: added_points) {
            index++;
            Point current = added_point->get_coord();
            if (pt0 == current) {
                b0 = false;
                p0 = added_point;
            }
            if (pt1 == current) {
                b1 = false;
                p1 = added_point;
            }
            if (pt2 == current) {
                b2 = false;
                p2 = added_point;
            }
        }

        if (b0) {
            p0 = std::make_shared<Mpoint>(pt0, mesh.nvertices() + count);
            count++;
        }
        if (b1) {
            p1 = std::make_shared<Mpoint>(pt1, mesh.nvertices() + count);
            count++;
        }
        if (b2) {
            p2 = std::make_shared<Mpoint>(pt2, mesh.nvertices() + count);
            count++;
        }

        if (b0) {
            mesh.push_point(p0);
            added_points.push_back(p0);
        }
        if (b1) {
            mesh.push_point(p1);
            added_points.push_back(p1);
        }
        if (b2) {
            mesh.push_point((p2));
            added_points.push_back(p2);
        }

        std::vector<int> tr_nbours(4);
        Triangle t0(p2, p0, p1, tot_triangles);
        tr_nbours.push_back(tot_triangles);
        tot_triangles++;
        Triangle t1(p1, v0, p2, tot_triangles);
        tr_nbours.push_back(tot_triangles);
        tot_triangles++;
        Triangle t2(p0, v2, p1, tot_triangles);
        tr_nbours.push_back(tot_triangles);
        tot_triangles++;
        Triangle t3(p2, v1, p0, tot_triangles);
        tr_nbours.push_back(tot_triangles);
        tot_triangles++;

        old_tr_nbours.push_back(tr_nbours);

        mesh.push_triangle(t0);
        mesh.push_triangle(t1);
        mesh.push_triangle(t2);
        mesh.push_triangle(t3);
    }

    for (auto i = mesh.vbegin(); i != mesh.vend(); i++)
        (*i)->normalize();
}

Mesh make_mesh_from_icosa(int n) {

    Mesh ret;

    const double tau = 0.8506508084;
    const double one = 0.5257311121;
    //creates regular icosahedron ////////////
    std::shared_ptr<Mpoint> ZA = std::make_shared<Mpoint>(tau, one, 0, 0);
    std::shared_ptr<Mpoint> ZB = std::make_shared<Mpoint>(-tau, one, 0, 1);
    std::shared_ptr<Mpoint> ZC = std::make_shared<Mpoint>(-tau, -one, 0, 2);
    std::shared_ptr<Mpoint> ZD = std::make_shared<Mpoint>(tau, -one, 0, 3);
    std::shared_ptr<Mpoint> YA = std::make_shared<Mpoint>(one, 0, tau, 4);
    std::shared_ptr<Mpoint> YB = std::make_shared<Mpoint>(one, 0, -tau, 5);
    std::shared_ptr<Mpoint> YC = std::make_shared<Mpoint>(-one, 0, -tau, 6);
    std::shared_ptr<Mpoint> YD = std::make_shared<Mpoint>(-one, 0, tau, 7);
    std::shared_ptr<Mpoint> XA = std::make_shared<Mpoint>(0, tau, one, 8);
    std::shared_ptr<Mpoint> XB = std::make_shared<Mpoint>(0, -tau, one, 9);
    std::shared_ptr<Mpoint> XC = std::make_shared<Mpoint>(0, -tau, -one, 10);
    std::shared_ptr<Mpoint> XD = std::make_shared<Mpoint>(0, tau, -one, 11);

    Triangle t0(YD, XA, YA, 0);
    Triangle t1(XB, YD, YA, 1);
    Triangle t2(XD, YC, YB, 2);
    Triangle t3(YC, XC, YB, 3);
    Triangle t4(ZD, YA, ZA, 4);
    Triangle t5(YB, ZD, ZA, 5);
    Triangle t6(ZB, YD, ZC, 6);
    Triangle t7(YC, ZB, ZC, 7);
    Triangle t8(XD, ZA, XA, 8);
    Triangle t9(ZB, XD, XA, 9);
    Triangle t10(ZD, XC, XB, 10);
    Triangle t11(XC, ZC, XB, 11);
    Triangle t12(ZA, YA, XA, 12);
    Triangle t13(YB, ZA, XD, 13);
    Triangle t14(ZD, XB, YA, 14);
    Triangle t15(XC, ZD, YB, 15);
    Triangle t16(ZB, XA, YD, 16);
    Triangle t17(XD, ZB, YC, 17);
    Triangle t18(XB, ZC, YD, 18);
    Triangle t19(ZC, XC, YC, 19);

    ret.push_point(ZA);
    ret.push_point(ZB);
    ret.push_point(ZC);
    ret.push_point(ZD);
    ret.push_point(YA);
    ret.push_point(YB);
    ret.push_point(YC);
    ret.push_point(YD);
    ret.push_point(XA);
    ret.push_point(XB);
    ret.push_point(XC);
    ret.push_point(XD);

    ret.push_triangle(t0);
    ret.push_triangle(t1);
    ret.push_triangle(t2);
    ret.push_triangle(t3);
    ret.push_triangle(t4);
    ret.push_triangle(t5);
    ret.push_triangle(t6);
    ret.push_triangle(t7);
    ret.push_triangle(t8);
    ret.push_triangle(t9);
    ret.push_triangle(t10);
    ret.push_triangle(t11);
    ret.push_triangle(t12);
    ret.push_triangle(t13);
    ret.push_triangle(t14);
    ret.push_triangle(t15);
    ret.push_triangle(t16);
    ret.push_triangle(t17);
    ret.push_triangle(t18);
    ret.push_triangle(t19);

    for (auto i = ret.tbegin(); i != ret.tend(); i++)
        i->swap_orientation();  //changes triangle orientation

    for (int io = 0; io < n; io++)
        retessellate(ret);

    std::vector<float> tmp_pvalues(ret.nvertices(), 0);
    ret.push_pvalues(tmp_pvalues);

    std::vector<float> tmp_tvalues(ret.ntriangles(), 0);
    ret.push_tvalues(tmp_tvalues);

    return ret;
}

void check_scale(Mesh& in, const Mesh& ref) {

    Point cr = in.get_coord(0),
          cr2 = in.get_coord(1),
          cr3 = ref.get_coord(1);

    if (std::abs(cr.norm() - cr2.norm()) > 1e-3 ||
        std::abs(cr.norm() - cr3.norm()) > 1e-3 ||
        std::abs(cr2.norm() - cr3.norm()) > 1e-3)
        true_rescale(in, cr3.norm());
}

void true_rescale(Mesh& mesh, double rad) {
// rescales sphere vertices to have equal radii
    for (int i = 0; i < mesh.nvertices(); i++)
    {
        Point cr = mesh.get_coord(i);
        cr.normalize();
        cr = cr * rad;
        mesh.set_coord(i, cr);
    }
}

void recentre(Mesh& sphere) {

    Point mean = sphere.estimate_origin();

    if(mean.norm() > 1e-2)
    {
        NEWMAT::Matrix translate(4, 4);
        translate = 0;
        translate(1, 1) = 1;
        translate(2, 2) = 1;
        translate(3, 3) = 1;
        translate(4, 4) = 1;
        translate(1, 4) = -mean.X;
        translate(2, 4) = -mean.Y;
        translate(3, 4) = -mean.Z;

        for (auto i = sphere.vbegin(); i != sphere.vend(); i++)
        {
            NEWMAT::ColumnVector P_in(4);
            Point p = (*i)->get_coord(), p2;
            if(p.norm() != 0.0)
            {
                P_in(1) = p.X;
                P_in(2) = p.Y;
                P_in(3) = p.Z;
                P_in(4) = 1;
                P_in = translate * P_in;
                p2.X = P_in(1);
                p2.Y = P_in(2);
                p2.Z = P_in(3);
                (*i)->set_coord(p2);
            }
        }
    }
}

bool check_for_intersections(int ind, double eps, Mesh& IN) {

    Point c = IN.get_coord(ind);
    Point V = c;
    V.normalize();
    int a = 0;

    const Triangle& tr = IN.get_triangle_from_vertex(ind, 0);
    c = IN.get_coord(ind);

    Point N = tr.normal();
    for (auto j = IN.tIDbegin(ind); j != IN.tIDend(ind); j++)
    {
        const Triangle& tr2 = IN.get_triangle(*j);
        Point N2 = tr2.normal();

        a = a || ((N|N2) <= 0.5);

        if (a == 1) break;
    }

    return a;
}

Point spatialgradient(int index, const Mesh& SOURCE) {
// spatial gradient for vertex trianglular mesh
    Point grad, ci = SOURCE.get_coord(index);

    for (auto j = SOURCE.tIDbegin(index); j != SOURCE.tIDend(index); j++)
    {
        Point v0 = SOURCE.get_triangle_vertex(*j, 0),
              v1 = SOURCE.get_triangle_vertex(*j, 1),
              v2 = SOURCE.get_triangle_vertex(*j, 2),
              dA;

        if ((ci - v0).norm() == 0)
            dA = computeGradientOfBarycentricTriangle(v1, v2, v0);
        else if ((ci - v1).norm() == 0)
            dA = computeGradientOfBarycentricTriangle(v2, v0, v1);
        else
            dA = computeGradientOfBarycentricTriangle(v0, v1, v2);

        grad = grad + dA;
    }
    return grad;
}

void unfold(Mesh& SOURCE) {

    // check for intersections
    std::vector<int> foldedvertices;
    std::vector<Point> foldinggradients;
    Point grad, ppp, ci;
    double current_stepsize;
    Mesh tmp = SOURCE;
    Point pp;
    bool folded = true;
    int it = 0;
    const double RAD = 100;

    while (folded)
    {
        foldedvertices.clear();
        foldinggradients.clear();

        for (int i = 0; i < SOURCE.nvertices(); i++)
            if (check_for_intersections(i, 0, SOURCE))
                // identfy points that are flipped.
                foldedvertices.push_back(i);

        if (foldedvertices.empty()) folded = false;
        else if (it % 100 == 0)
            cout << " mesh is folded, total folded vertices =" << foldedvertices.size() << " unfold " << endl;

        for (int foldedvertice : foldedvertices)
        {
            grad = spatialgradient(foldedvertice, SOURCE);  // estimates gradient of triangle area for all vertices that are detected as folded
            foldinggradients.push_back(grad);
        }

        for (unsigned int i = 0; i < foldedvertices.size(); i++)
        {
            current_stepsize = 1;
            ci = SOURCE.get_coord(foldedvertices[i]);
            grad = foldinggradients[i];
            do
            {
                pp = ci - grad * current_stepsize; // gradient points in direction of increasing area therefore opposite to desired direction
                pp.normalize();
                SOURCE.set_coord(foldedvertices[i],pp*RAD);
                current_stepsize /= 2.0;
            }
            while(check_for_intersections(foldedvertices[i],0,SOURCE) == 1 && current_stepsize > 1e-3);

            SOURCE.set_coord(foldedvertices[i],pp*RAD);
        }
        it++;
        if(it == 1000) break;
    }
}

NEWMAT::Matrix get_coordinate_transformation(double dNdT1,double dNdT2, NEWMAT::ColumnVector& Norm) {

    Point G1(1, 0, dNdT1);
    Point G2(0, 1, dNdT2);
    Point G3 = G1 * G2;

    G3 = G3 / sqrt((G3|G3));

    NEWMAT::Matrix G(3,3);

    G(1,1) = G1.X;
    G(1,2) = G2.X;
    G(1,3) = G3.X;
    G(2,1) = G1.Y;
    G(2,2) = G2.Y;
    G(2,3) = G3.Y;
    G(3,1) = G1.Z;
    G(3,2) = G2.Z;
    G(3,3) = G3.Z;
    Norm(1) = G3.X;
    Norm(2) = G3.Y;
    Norm(3) = G3.Z;

    return (G.i()).t();
}

NEWMAT::ColumnVector calculate_strains(int index, const std::vector<int>& kept, const Mesh& orig, const Mesh& final, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches){

    NEWMAT::ColumnVector STRAINS(4); STRAINS = 0;
    Point Normal_O = orig.get_normal(index - 1);
    NEWMAT::Matrix TRANS(3,3),
                   principal_STRAINS(2, orig.nvertices());
    double T1_coord = 0.0, T2_coord = 0.0;
    Point orig_tang,final_tang,tmp;
    //Tangent Tang;
    NEWMAT::Matrix a, b, c, d;
    int maxind, minind;

    Tangs T = calculate_tangs(index - 1, orig);

    TRANS = form_matrix_from_points(T.e1, T.e2, Normal_O, true);

    NEWMAT::Matrix pseudo_inv, pseudo_V, pseudo_U;
    NEWMAT::DiagonalMatrix pseudo_D;
    NEWMAT::Matrix alpha(kept.size(),6),
                   N(kept.size(),1),
                   n(kept.size(),1),
                   t1(kept.size(),1),
                   t2(kept.size(),1);
    alpha = 0; N = 0; n = 0; t1 = 0; t2 = 0;

    for(int j = 0; j < (int)kept.size(); j++)
    {
        tmp = orig.get_coord(kept[j]) - orig.get_coord(index - 1);
        projectPoint(tmp, T, T1_coord, T2_coord);
        N(j+1,1) = tmp | Normal_O;

        // fit poly
        alpha(j+1,2) = T1_coord;
        alpha(j+1,3) = T2_coord;
        alpha(j+1,4) = 0.5 * T1_coord * T1_coord;
        alpha(j+1,5) = 0.5 * T2_coord * T2_coord;
        alpha(j+1,6) = T1_coord * T2_coord;

        // get transformed coords
        tmp = final.get_coord(kept[j]) - final.get_coord(index - 1);
        projectPoint(tmp, T, T1_coord, T2_coord);

        n(j + 1,1) = tmp | Normal_O;
        t1(j + 1,1) = T1_coord;
        t2(j + 1,1) = T2_coord;
    }

    SVD(alpha, pseudo_D, pseudo_U, pseudo_V);

    for (int i = 1; i <= pseudo_D.Nrows(); i++)
        if (pseudo_D(i) != 0)
            pseudo_D(i) = 1 / pseudo_D(i);

    pseudo_inv = pseudo_V * pseudo_D * pseudo_U.t();

    // get coeffieicnts for different fits
    a = pseudo_inv * N;
    b = pseudo_inv * t1;
    c = pseudo_inv * t2;
    d = pseudo_inv * n;

    double dNdT1, dNdT2, dt1dT1, dt1dT2, dt2dT1, dt2dT2, dndT1, dndT2;
    dNdT1 = a(2,1); dNdT2 = a(3,1);
    dt1dT1 = b(2,1); dt1dT2 = b(3,1);
    dt2dT1 = c(2,1); dt2dT2 = c(3,1);
    dndT1 = d(2,1); dndT2 = d(3,1);

    NEWMAT::ColumnVector G3(3);
    NEWMAT::Matrix G_cont = get_coordinate_transformation(dNdT1, dNdT2, G3);

    Point g1(dt1dT1,dt2dT1,dndT1);
    Point g2(dt1dT2,dt2dT2,dndT2);
    Point g3 = g1 * g2; g3 = g3 / sqrt((g3|g3));

    NEWMAT::Matrix g(3,3);
    g(1,1)=g1.X; g(1,2)=g2.X; g(1,3)=g3.X;
    g(2,1)=g1.Y; g(2,2)=g2.Y; g(2,3)=g3.Y;
    g(3,1)=g1.Z; g(3,2)=g2.Z; g(3,3)=g3.Z;

    NEWMAT::Matrix F = g * G_cont.t();
    NEWMAT::Matrix C = F.t() * F;

    NEWMAT::DiagonalMatrix Omega;
    NEWMAT::Matrix U, Umax, Umin;
    SVD(C,Omega,U);

    NEWMAT::Matrix mm = G3.t() * U;

    if (abs(mm(1, 1)) >= abs(mm(1, 2)) && abs(mm(1, 1)) >= abs(mm(1, 3)))
        if (sqrt(Omega(2)) > sqrt(Omega(3)))
        {
            maxind = 2;
            minind = 3;
        }
        else
        {
            maxind = 3;
            minind = 2;
        }

    else if (abs(mm(1, 2)) >= abs(mm(1, 1)) && abs(mm(1, 2)) >= abs(mm(1, 3)))
        if (sqrt(Omega(1)) > sqrt(Omega(3)))
        {
            maxind = 1;
            minind = 3;
        }
        else
        {
            maxind = 3;
            minind = 1;
        }
    else
    {
        if (sqrt(Omega(1)) > sqrt(Omega(2)))
        {
            maxind = 1;
            minind = 2;
        }
        else
        {
            maxind = 2;
            minind = 1;
        }
    }

    STRAINS(1) = sqrt(Omega(maxind));
    STRAINS(2) = sqrt(Omega(minind));
    Umax = TRANS.i() * U.SubMatrix(1,3,maxind,maxind);
    Umin = TRANS.i() * U.SubMatrix(1,3,minind,minind);

    if(PrincipalStretches)
    {
        (*PrincipalStretches)(index,1) = Umax(1,1);
        (*PrincipalStretches)(index,2) = Umax(2,1);
        (*PrincipalStretches)(index,3) = Umax(3,1);
        (*PrincipalStretches)(index,4) = Umin(1,1);
        (*PrincipalStretches)(index,5) = Umin(2,1);
        (*PrincipalStretches)(index,6) = Umin(3,1);
    }

    STRAINS(3) = 0.5 * (STRAINS(1)*STRAINS(1) - 1);
    STRAINS(4) = 0.5 * (STRAINS(2)*STRAINS(2) - 1);

    return STRAINS;
}

Mesh calculate_strains(double fit_radius, const Mesh& orig, const Mesh &final, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches){

    Mesh strain = final;
    NEWMAT::Matrix strains(4, orig.nvertices());
    NEWMAT::ColumnVector node_strain;
    double dir_chk1, dir_chk2, dist, fit_temp;

    vector<int> kept;

    int neighbour, numneighbours = orig.nvertices();

    if(PrincipalStretches) PrincipalStretches->ReSize(orig.nvertices(), 6);

    for (int index = 1; index <= orig.nvertices(); index++)
    {
        //if(REL) numneighbours = REL->Nrows(index);
        fit_temp = fit_radius;
        while (kept.size() <= 8)
        {
            kept.clear();
            for(int j = 1; j <= numneighbours; j++)
            {
                /*if(REL) neighbour = (*REL)(j,index);
                else*/ neighbour = j;

                dist = (orig.get_coord(index - 1) - orig.get_coord(neighbour - 1)).norm();

                // reject points with normals in opposite directions
                dir_chk1 = orig.get_normal(neighbour - 1) | orig.get_normal(index - 1);
                dir_chk2 = 1; //Normals_F[j-1]|Normals_F[index-1]; ////SET to 1?

                if(dist <= fit_temp && dir_chk1 >= 0 && dir_chk2 >= 0)
                    kept.push_back(neighbour - 1);
            }

            if(kept.size() > 8)
            {
                node_strain = calculate_strains(index, kept, orig, final, PrincipalStretches);
                strains(1, index) = node_strain(1);
                strains(2, index) = node_strain(2);
                strains(3, index) = node_strain(3);
                strains(4, index) = node_strain(4);
            }
            else
                fit_temp += 0.5;
        }
        kept.clear();
    }

    strain.set_pvalues(strains);

    return strain;
}

Tangs calculate_tangs(int ind, const Mesh& SPH_in) {

    Tangs T;
    double mag;
    Point a = SPH_in.local_normal(ind);
    Point tmp = SPH_in.get_coord(ind);

    if((a|tmp) < 0) a = a * -1;

    if(abs(a.X) >= abs(a.Y) && abs(a.X) >= abs(a.Z))
    {
        mag = sqrt(a.Z*a.Z + a.Y*a.Y);
        if(mag == 0)
        {
            T.e1.X = 0;
            T.e1.Y = 0;
            T.e1.Z = 1;
        }
        else
        {
            T.e1.X = 0;
            T.e1.Y = -a.Z /mag;
            T.e1.Z = a.Y /mag;
        }
    }
    else if(abs(a.Y) >= abs(a.X) && abs(a.Y) >= abs(a.Z))
    {
        mag = sqrt(a.Z*a.Z + a.X*a.X);
        if(mag == 0)
        {
            T.e1.X = 0;
            T.e1.Y = 0;
            T.e1.Z = 1;
        }
        else
        {
            T.e1.X = -a.Z/mag;
            T.e1.Y = 0;
            T.e1.Z = a.X /mag;
        }
    }
    else
    {
        mag = sqrt(a.Y*a.Y + a.X*a.X);
        if(mag == 0)
        {
            T.e1.X = 1;
            T.e1.Y = 0;
            T.e1.Z = 0;
        }
        else
        {
            T.e1.X = -a.Y/mag;
            T.e1.Y = a.X/mag;
            T.e1.Z = 0;
        }
    }
    T.e2 = a * T.e1;
    T.e2.normalize();

    return T;
}

double compute_vertex_area(int ind, const Mesh& mesh) {

    double sum = 0;

    for (auto i = mesh.tIDbegin(ind); i != mesh.tIDend(ind); i++)
        sum += mesh.get_triangle_area(*i);

    return sum / mesh.get_total_triangles(ind);
}

NEWMAT::ReturnMatrix rotate_vec(const NEWMAT::ColumnVector& vec, double w1, double w2, double w3) {

    NEWMAT::Matrix rotation(3, 3);
    NEWMAT::ColumnVector rotated_vec(3);

    rotation << cos(w2) * cos(w3)
        << -cos(w1) * sin(w3) + sin(w1) * sin(w2) * cos(w3)
        << sin(w1) * sin(w3) + cos(w1) * sin(w2) * cos(w3)
        << cos(w2) * sin(w3)
        << cos(w1) * cos(w3) + sin(w1) * sin(w2) * sin(w3)
        << -sin(w1) * cos(w3) + cos(w1) * sin(w2)*sin(w3)
        << -sin(w2)
        << sin(w1) * cos(w2)
        << cos(w1) * cos(w2);

    rotation = rotation.t();

    rotated_vec = rotation * vec;
    rotated_vec.Release();

    return rotated_vec;
}

bool operator==(const Mesh& M1, const Mesh& M2) {
    if (M1.nvertices() != M2.nvertices()) return false;
    for (int i = 0; i < M1.nvertices(); i++)
        if (M1.get_point(i) != M2.get_point(i)) return false;
    return true;
}

} //namespace newresampler
