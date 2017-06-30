/* This source file contains the Palander simulation engine, to be used in 
 * conjunction with the Blender add-on palander_addon.py and the scripts
 * palander_remoterun.sh and palander_scriptbase. All files are available at
 * <https://github.com/FilmakademieRnD/Palander>
 *
 * Copyright (c) 2017 Animationsinstitut, Filmakademie Baden-Wuerttemberg
 * <http://research.animationsinstitut.de/>
 *
 * This project has been realized in the scope of the Media Solution Center
 * project, funded by the state of Baden-Wuerttemberg, Germany:
 * <http://msc-bw.de/en>
 *
 * This file makes use of the Palabos library, which is free software released
 * under the GNU Affero General Public License (AGPL) as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version. Therefore, this file is also released under the GNU AGPL;
 * you can redistribute it and/or modify it under the terms of the GNU AGPL.
 *
 * The most recent official release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * DISCLAIMER: Palander and the library Palabos are free software.
 *
 * Palander and the library are distributed in the hope that they will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero
 * General Public License for more details.
 *
 * You may have received a copy of the GNU Affero General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos3D.h"
#include "palabos3D.hh"

#include <zlib.h>

using namespace plb;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

#define OUTDIR     "tmp/"
#define INFILE     "toPalabos.xml"

#define FLAG_SAVE  "tmp/flag"
#define J_SAVE     "tmp/j"
#define LAT_SAVE   "tmp/lattice"
#define MASS_SAVE  "tmp/mass"
#define NORM_SAVE  "tmp/normal"
#define RHO_SAVE   "tmp/rhoBar"
#define VF_SAVE    "tmp/volFrac"

#define BDNUM      1000
#define CTANG      80.0

#define FPS        25
#define PADDING    5
#define STRMAX     500
#define MAX_REFINE 16

typedef double T;


// Smagorinsky constant for LES model.
const T cSmago = 0.14;

const T rhoEmpty = T(1);
    
plint iT, maxIter, outIter, statIter, saveOffset, saveIndex;
plint objCount, fluidCount, staticSurfaces, movingSurfaces;
plint N, nx, ny, nz;
T lx, ly, lz, size, lxHalf, lyHalf, lzHalf;
T dx, dy, dz, dxmin, dymin, dzmin, dxmax, dymax, dzmax;
T delta_t, delta_x, rotNorm, velNorm;
T uLB, nuPhys, nuLB, tau, omega, Bo, surfaceTensionLB, contactAngle;
bool disableST, continuous, smooooth, useMask, foundMask, toCluster;
bool BOBJpre, BOBJfin, STLint, STLsmint, STLmov, VTKvf, VTKvel, VTKvort;
Array<T,3> externalForce;
std::string basePath, cache, maskFile, *objType;
std::vector<float> location, rotation, sceneScale;
std::vector<T> xmin, ymin, zmin, xmax, ymax, zmax;
std::vector<T> *matrix, *vertices;
std::vector<T> moveX, moveY, moveZ, rotX, rotY, rotZ;
Array<T,3> moveStep, rotStep;

std::string outDir = OUTDIR;

std::vector<DEFscaledMesh<T>* > fluidMesh;
std::vector<DEFscaledMesh<T>* > solidMesh;

MultiScalarField3D<int> *mask;
int ***outMask;

void setupParameters() {
    delta_x = size / N;
    delta_t = delta_x * uLB;
    nx = util::roundToInt(lx / delta_x);
    ny = util::roundToInt(ly / delta_x);
    nz = util::roundToInt(lz / delta_x);

    // Gravity in lattice units.
    T gLB = 9.8 * delta_t * delta_t/delta_x;
    externalForce = Array<T,3>(0., 0., -gLB);
    tau            = (nuPhys*DESCRIPTOR<T>::invCs2*delta_t)/(delta_x*delta_x) + 0.5;
    omega          = 1./tau;    
    nuLB           = (tau-0.5)*DESCRIPTOR<T>::cs2; // Viscosity in lattice units.
    
    surfaceTensionLB = rhoEmpty * gLB * N * N / Bo;

    lxHalf = lx/2.0; lyHalf = ly/2.0; lzHalf = lz/2.0;

    if (useMask) mask = new MultiScalarField3D<int>(nx, ny, nz);

}

int initialFluidFlags(plint x, plint y, plint z)
{
    plint i, index, fluidIntersections, solidIntersections, checkIntersections, numTriangles;

    // Rules for origin location:
    //  - MUST BE outside of simulation cell (at least one coordinate negative or greater than domain size)
    //  - AVOID GLITCHES: clear view to each point in cell (only one coordinate negative or greater than domain size)

    Array<T,3> origin = Array<T,3>(lxHalf,lyHalf,-lzHalf);
    Array<T,3> checkOrigin = Array<T,3>(lxHalf,lyHalf,lz+lzHalf);
    Array<T,3> backupOrigin = Array<T,3>(-lxHalf,lyHalf,lzHalf);
    Array<T,3> point = Array<T,3>(x*delta_x,y*delta_x,z*delta_x);

    fluidIntersections = solidIntersections = checkIntersections = 0;
    for (i = 0; i < fluidCount; i++) {
        TriangularSurfaceMesh<T> tempMesh = fluidMesh[i]->getMesh();
        numTriangles = tempMesh.getNumTriangles();
        for (index = 0; index < numTriangles; index++) {
            fluidIntersections += tempMesh.segmentIntersectsTriangle(point, origin, index) ? 1 : 0;
            checkIntersections += tempMesh.segmentIntersectsTriangle(point, checkOrigin, index) ? 1 : 0;
        }
        if (fluidIntersections != checkIntersections) { // Necessary check to remove bug
#ifdef PLB_DEBUG
            pcout << "WARNING: fluid object may have been read incorrectly - please check the results!" << std::endl;
#endif
            for (index = 0, fluidIntersections = 0; index < numTriangles; index++)
                fluidIntersections += tempMesh.segmentIntersectsTriangle(point, backupOrigin, index) ? 1 : 0;
        }
        if (fluidIntersections % 2) break;
    }
    for (i = 0; i < staticSurfaces; i++) {
        TriangularSurfaceMesh<T> tempMesh = solidMesh[i]->getMesh();
        numTriangles = tempMesh.getNumTriangles();
        for (index = 0; index < numTriangles; index++) {
            solidIntersections += tempMesh.segmentIntersectsTriangle(point, origin, index) ? 1 : 0;
            checkIntersections += tempMesh.segmentIntersectsTriangle(point, checkOrigin, index) ? 1 : 0;
        }
        if (solidIntersections != checkIntersections) { // Necessary check to remove bug
#ifdef PLB_DEBUG
            pcout << "WARNING: solid object may have been read incorrectly - please check the results!" << std::endl;
#endif
            for (index = 0, solidIntersections = 0; index < numTriangles; index++)
                solidIntersections += tempMesh.segmentIntersectsTriangle(point, backupOrigin, index) ? 1 : 0;
        }
        if (solidIntersections % 2) break;
    }

    if (solidIntersections % 2) {
        if (useMask) outMask[x][y][z] = freeSurfaceFlag::wall;
        return freeSurfaceFlag::wall;
    } else if (fluidIntersections % 2) {
        if (useMask) outMask[x][y][z] = freeSurfaceFlag::fluid;
        return freeSurfaceFlag::fluid;
    } else {
        if (useMask) outMask[x][y][z] = freeSurfaceFlag::empty;
        return freeSurfaceFlag::empty;
    }

}

Array<T,3> calculateVelocity(plint frame, Array<T,3> vertex, const std::vector<T> &posX, const std::vector<T> &posY, const std::vector<T> &posZ, const std::vector<T> &rotX, const std::vector<T> &rotY, const std::vector<T> &rotZ) {
    Array<T,3> oldPosition = Array<T,3>(posX[frame],posY[frame],posZ[frame]);
    Array<T,3> newPosition = Array<T,3>(posX[frame+1],posY[frame+1],posZ[frame+1]);
    Array<T,3> oldRotation = Array<T,3>(rotX[frame],rotY[frame],rotZ[frame]);
    Array<T,3> newRotation = Array<T,3>(rotX[frame+1],rotY[frame+1],rotZ[frame+1]);
    Array<T,3> pointOnAxis = Array<T,3>(nx*(posX[frame]-dxmin)/dx,ny*(posY[frame]-dymin)/dy,nz*(posZ[frame]-dzmin)/dz);
    Array<T,3> totalStep = Array<T,3>(0.0,0.0,0.0);

    if (frame < rotX.size()-1) {
        Array<T,3> angularVelocity = rotNorm * (newRotation - oldRotation);
        if (angularVelocity != Array<T,3>(0.0,0.0,0.0)) {
            Array<T,3> rotationAxis = angularVelocity / norm<T,3>(angularVelocity);
            totalStep += getRotationalVelocity(vertex, angularVelocity, rotationAxis, pointOnAxis);
        }
    }

    if (frame < posX.size()-1) {
        totalStep += velNorm * (newPosition - oldPosition);
    }

    return totalStep;

}

bool checkInlet(plint frame, std::vector<plint> enabledFrame, std::vector<plint> enabledSwitch)
{
    plint i, iMax = enabledFrame.size();

    if (iMax == 0) return(false);
    for (i = 0; i < iMax && frame >= enabledFrame[i]; i++);
    if (i > 0) --i;

    return(enabledSwitch[i] == 1 ? true : false);

}

void writeResults(FreeSurfaceFields3D<T,DESCRIPTOR> *fields, plint iT)
{
    std::auto_ptr<MultiScalarField3D<T> > smoothVF(lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction,
                fields->volumeFraction.getBoundingBox().enlarge(-1)));
    std::auto_ptr<MultiScalarField3D<T> > smoooothVF(lbmSmoothen<T,DESCRIPTOR>(*lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction, fields->volumeFraction.getBoundingBox().enlarge(-1)),
                fields->volumeFraction.getBoundingBox()));

    std::vector<T> isoLevels;
    isoLevels.push_back(0.5);
    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;

    // Preview images
    if (smooooth)
        isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox());
    else
        isoSurfaceMarchingCube(triangles, fields->volumeFraction, isoLevels, fields->volumeFraction.getBoundingBox());
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(delta_x);
        // STL output filename indexed by milliseconds in simulation time
        if (STLint) triangleSet.writeBinarySTL(createFileName(outDir + "interface_", (int)(iT*1000.0*delta_t), PADDING)+"ms.stl");
        if (BOBJpre) {
            // The following scale and translate are necessary for Blender coordinate conversion
            triangleSet.scale(dx/lx,dy/ly,dz/lz);
            triangleSet.translate(Array<T,3>(dxmin,dymin,dzmin));

            // Blender bobj.gz output
            if (global::mpi().isMainProcessor() && saveIndex >= 0) {
                int numVertices, tempI;
                float tempF;
                std::string outfile;
                gzFile gzf;
                if (toCluster)
                    outfile = createFileName(outDir + "fluidsurface_preview_", saveIndex, 4) + ".bobj.gz";
                else
                    outfile = createFileName(cache + "/" + "fluidsurface_preview_", saveIndex, 4) + ".bobj.gz";
                DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(triangleSet);
                std::vector<Array<T,3> > vertexList = defMesh->getVertexList();
                std::vector<Edge> edgeList = defMesh->getEdgeList();
                delete defMesh;

                numVertices = (int)(vertexList.size());
                gzf = gzopen(outfile.c_str(), "wb1");

                // Write vertex list
                gzwrite(gzf, &numVertices, sizeof(int));
                for (size_t i = 0; i < numVertices; i++) {
                    for (int j = 0; j < 3; j++) {
                        tempF = ((float)(vertexList[i][j])-location[j])/sceneScale[j];
                        gzwrite(gzf, &tempF, sizeof(float));
                    }
                }
                // Write (dummy) normals list
                tempF = 1.0;
                gzwrite(gzf, &numVertices, sizeof(int));
                for (size_t i = 0; i < numVertices; i++) {
                    for (int j = 0; j < 3; j++) {
                        gzwrite(gzf, &tempF, sizeof(float));
                    }
                }
                // Write triangle index list
                tempI = triangles.size();
                gzwrite(gzf, &tempI, sizeof(int));
                for (plint i = 0; i < tempI; i++) {
                    for (plint j = 2; j >= 0; j--) {
                        int index = edgeList[3*i + j].pv;
                        gzwrite(gzf, &index, sizeof(int));
                    }
                }
                gzclose(gzf);
            }
        }
    }

    // Smoothened final images
    triangles.clear();
    if (smooooth)
        isoSurfaceMarchingCube(triangles, *smoooothVF, isoLevels, smoooothVF->getBoundingBox());
    else
        isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox());
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(delta_x);
        // STL output filename indexed by milliseconds in simulation time
        if (STLsmint) triangleSet.writeBinarySTL(createFileName(outDir + "smoothedInterface_", (int)(iT*1000.0*delta_t), PADDING)+"ms.stl");
        if (BOBJfin) {
            // The following scale and translate are necessary for Blender coordinate conversion
            triangleSet.scale(dx/lx,dy/ly,dz/lz);
            triangleSet.translate(Array<T,3>(dxmin,dymin,dzmin));

            // Blender bobj.gz output
            if (global::mpi().isMainProcessor() && saveIndex >= 0) {
                int numVertices, tempI;
                float tempF;
                std::string outfile;
                gzFile gzf;
                if (toCluster)
                    outfile = createFileName(outDir + "fluidsurface_final_", saveIndex, 4) + ".bobj.gz";
                else
                    outfile = createFileName(cache + "/" + "fluidsurface_final_", saveIndex, 4) + ".bobj.gz";
#ifdef PLB_DEBUG
                pcout << "Saving to file " << outfile << std::endl;
#endif
                DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(triangleSet);
                std::vector<Array<T,3> > vertexList = defMesh->getVertexList();
                std::vector<Edge> edgeList = defMesh->getEdgeList();
                delete defMesh;

                numVertices = (int)(vertexList.size());
                gzf = gzopen(outfile.c_str(), "wb1");

                // Write vertex list
                gzwrite(gzf, &numVertices, sizeof(int));
                for (size_t i = 0; i < numVertices; i++) {
                    for (int j = 0; j < 3; j++) {
                        tempF = ((float)(vertexList[i][j])-location[j])/sceneScale[j];
                        gzwrite(gzf, &tempF, sizeof(float));
                    }
                }
                // Write (dummy) normals list
                tempF = 1.0;
                gzwrite(gzf, &numVertices, sizeof(int));
                for (size_t i = 0; i < numVertices; i++) {
                    for (int j = 0; j < 3; j++) {
                        gzwrite(gzf, &tempF, sizeof(float));
                    }
                }
                // Write triangle index list
                tempI = triangles.size();
                gzwrite(gzf, &tempI, sizeof(int));
                for (plint i = 0; i < tempI; i++) {
                    for (plint j = 2; j >= 0; j--) {
                        int index = edgeList[3*i + j].pv;
                        gzwrite(gzf, &index, sizeof(int));
                    }
                }
                gzclose(gzf);
            }
        }
    }

    if (VTKvf || VTKvel || VTKvort) {
        VtkImageOutput3D<T> vtkOut(createFileName("VTKoutput_", (int)(iT*1000.0*delta_t), PADDING)+"ms", 1.);
        if (VTKvf) vtkOut.writeData<float>(fields->volumeFraction, "vf", 1.);
        if (VTKvel) vtkOut.writeData<3,float>(*computeVelocity(fields->lattice), "vel", delta_x/delta_t);
        if (VTKvort) vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(fields->lattice)), "vort", 1./delta_t);
    }

}

void getDimensions(plint index, std::vector<T> *matrix, std::vector<T> *vertices) {
    T xmult, xoff, ymult, yoff, zmult, zoff, x, y, z;

    for (plint k = 0; k < objCount; k++) {
        // Read critical scene object coordinates...
        xmult = matrix[k][0]; xoff = matrix[k][3];
        ymult = matrix[k][5]; yoff = matrix[k][7];
        zmult = matrix[k][10]; zoff = matrix[k][11];

        // ... and use them to find the bounding box of each object
        xmin.push_back(vertices[k][0]*xmult + xoff);
        xmax.push_back(vertices[k][0]*xmult + xoff);
        ymin.push_back(vertices[k][1]*ymult + yoff);
        ymax.push_back(vertices[k][1]*ymult + yoff);
        zmin.push_back(vertices[k][2]*zmult + zoff);
        zmax.push_back(vertices[k][2]*zmult + zoff);
        for (plint i = 3; i < vertices[k].size(); i+=3) {
            x = vertices[k][i]*xmult + xoff;
            y = vertices[k][i+1]*ymult + yoff;
            z = vertices[k][i+2]*zmult + zoff;
            if(x < xmin[k]) xmin[k] = x;
            else if(x > xmax[k]) xmax[k] = x;
            if(y < ymin[k]) ymin[k] = y;
            else if(y > ymax[k]) ymax[k] = y;
            if(z < zmin[k]) zmin[k] = z;
            else if(z > zmax[k]) zmax[k] = z;
        }
    }

}

class ForwardVelocity {
public:
    ForwardVelocity(std::vector<Array<T,3> > const& vertices_)
        : vertices(vertices_)
    { }
    Array<T,3> operator()(pluint id)
    {
        return moveStep;
    }
private:
    std::vector<Array<T,3> > const& vertices;

};

int main(int argc, char **argv)
{
    plbInit(&argc, &argv);
    global::directories().setInputDir("./");
    global::directories().setOutputDir(outDir+"/");

    plint domainIndex, inletIndex, outletIndex;
    const plint ibIter = 4;
    plint iniIter, iniTime, maxTime;
    std::string obj, dat, *fName;
    std::stringstream sstm;
    bool isInlet = false, isOutlet = false;
    std::vector<Box3D> bounds;
    std::vector<T> inVel;
    std::vector<Array<T,3> > velocity;
    std::vector<plint> enabledFrame, enabledSwitch;

    // Immersed Boundary Method variables
    std::vector<plint> startIds;
    std::vector<Array<T,3> > vertexArray;
    std::vector<T> areas;
    std::vector<int> flags;
    std::vector<ConnectedTriangleSet<T> > movingObstacles;

    Bo = BDNUM; contactAngle = CTANG;
    disableST = false; continuous = false;

#ifdef PLB_DEBUG
    pcout << "Reading Blender parameters..." << std::endl;
#endif
    XMLreader xmlFile(INFILE);
    xmlFile["Numerics"]["viscosity"].read(nuPhys);
    xmlFile["Numerics"]["maxTime"].read(maxTime);
    xmlFile["Numerics"]["resolution"].read(N);
    xmlFile["Numerics"]["saveOffset"].read(saveOffset);
    xmlFile["Numerics"]["latticeVel"].read(uLB);
    xmlFile["Numerics"]["disableST"].read(disableST);
    xmlFile["Numerics"]["continuous"].read(continuous);
    xmlFile["Numerics"]["extraSmooth"].read(smooooth);
    xmlFile["Numerics"]["useMask"].read(useMask);

    xmlFile["Output"]["maskFile"].read(maskFile);
    xmlFile["Output"]["basePath"].read(basePath);
    xmlFile["Output"]["toCluster"].read(toCluster);
    xmlFile["Output"]["BOBJpre"].read(BOBJpre);
    xmlFile["Output"]["BOBJfin"].read(BOBJfin);
    xmlFile["Output"]["STLint"].read(STLint);
    xmlFile["Output"]["STLsmint"].read(STLsmint);
    xmlFile["Output"]["STLmov"].read(STLmov);
    xmlFile["Output"]["VTKvf"].read(VTKvf);
    xmlFile["Output"]["VTKvel"].read(VTKvel);
    xmlFile["Output"]["VTKvort"].read(VTKvort);

    xmlFile["Geometry"]["objectCount"].read(objCount);
    bounds.resize(objCount);
    objType = new std::string[objCount];
    fName = new std::string[objCount];
    matrix = new std::vector<T>[objCount];
    vertices = new std::vector<T>[objCount];

    maxTime *= 1000.0;
    foundMask = false;
    fluidCount = staticSurfaces = movingSurfaces = 0;
    inVel.resize(0); velocity.resize(0); enabledFrame.resize(0); enabledSwitch.resize(0);
    fluidMesh.resize(0); solidMesh.resize(0);

    vertexArray.resize(0);
    areas.resize(0);
    flags.resize(0);
    startIds.resize(0);
    moveX.resize(0);
    moveY.resize(0);
    moveZ.resize(0);

    for (plint i = 0; i < objCount; i++) {
        Array<T,3> tempA;
        tempA.resetToZero();
        sstm << "object" << i;
        obj = sstm.str();
        sstm << ".stl";
        fName[i] = sstm.str();
        sstm.str(""); sstm.clear();
        sstm << "type" << i;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        xmlFile["Geometry"][obj][dat].read(objType[i]);
        sstm << "m" << i;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        xmlFile["Geometry"][obj][dat].read(matrix[i]);
        sstm << "v" << i;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        xmlFile["Geometry"][obj][dat].read(vertices[i]);
        if (objType[i] == "domain") {
            domainIndex = i;
            xmlFile["Geometry"][obj]["size"].read(size);
            xmlFile["Geometry"][obj]["location"].read(location);
            xmlFile["Geometry"][obj]["rotation"].read(rotation);
            xmlFile["Geometry"][obj]["scale"].read(sceneScale);
            xmlFile["Geometry"][obj]["cache"].read(cache);
            if (cache[0] == '/' && cache[1] == '/') {
                if (cache[2] == '/') cache.erase(0,3);
                else cache.erase(0,2);
                cache = basePath + cache;
            }
        }
        if (objType[i] == "fluid") {
            xmlFile["Geometry"][obj]["initialVelocity"].read(inVel);
            for (plint xyz = 0; xyz < 3; xyz++)
                tempA[xyz] = inVel[xyz];
            inVel.clear(); inVel.resize(0);
        }
        if (objType[i] == "obstacle") {
            try {
                xmlFile["Geometry"][obj]["motionPathX"].read(moveX);
                xmlFile["Geometry"][obj]["motionPathY"].read(moveY);
                xmlFile["Geometry"][obj]["motionPathZ"].read(moveZ);
                xmlFile["Geometry"][obj]["rotationPathX"].read(rotX);
                xmlFile["Geometry"][obj]["rotationPathY"].read(rotY);
                xmlFile["Geometry"][obj]["rotationPathZ"].read(rotZ);
            } catch (plb::PlbIOException) {
                objType[i] = "wall";
#ifdef PLB_DEBUG
                pcout << "Changed object " << i << " to type " << objType[i] << std::endl;
#endif
            }
        }
        if (objType[i] == "inflow") {
            isInlet = true;
            inletIndex = i;
            xmlFile["Geometry"][obj]["inflowVelocity"].read(inVel);
            for (plint xyz = 0; xyz < 3; xyz++)
                tempA[xyz] = inVel[xyz];
            inVel.clear(); inVel.resize(0);
            try {
                xmlFile["Geometry"][obj]["enabledFrame"].read(enabledFrame);
                xmlFile["Geometry"][obj]["enabledSwitch"].read(enabledSwitch);
            } catch (plb::PlbIOException) {
                enabledFrame.push_back(0);
                enabledSwitch.push_back(1);
            }
        }
        if (objType[i] == "outflow") {
            isOutlet = true;
            outletIndex = i;
        }
        velocity.push_back(tempA);
    }

    // Define Palabos simulation cell dimensions (lx, ly, lz)
    getDimensions(domainIndex, matrix, vertices);
    dxmin = xmin[domainIndex]; dymin = ymin[domainIndex]; dzmin = zmin[domainIndex];
    dxmax = xmax[domainIndex]; dymax = ymax[domainIndex]; dzmax = zmax[domainIndex];
    dx = dxmax-dxmin; dy = dymax-dymin; dz = dzmax-dzmin;
    if (dx >= dy) {
        if (dx >= dz) {
            lx = size;
            ly = size * dy/dx;
            lz = size * dz/dx;
        } else {
            lz = size;
            lx = size * dx/dz;
            ly = size * dy/dz;
        }
    } else {
        if (dy >= dz) {
            ly = size;
            lx = size * dx/dy;
            lz = size * dz/dy;
        } else {
            lz = size;
            lx = size * dx/dz;
            ly = size * dy/dz;
        }
    }

    // Convert min and max values into Palabos coordinates
    for (plint k = 0; k < objCount; k++) {
        xmin[k] = (xmin[k] - dxmin) * lx/dx;
        xmax[k] = (xmax[k] - dxmin) * lx/dx;
        ymin[k] = (ymin[k] - dymin) * ly/dy;
        ymax[k] = (ymax[k] - dymin) * ly/dy;
        zmin[k] = (zmin[k] - dzmin) * lz/dz;
        zmax[k] = (zmax[k] - dzmin) * lz/dz;
    }

#ifdef PLB_DEBUG
    pcout << "Setting up parameters..." << std::endl;
#endif
    setupParameters();
    outIter = util::roundToInt(1.0/(FPS*delta_t));
    statIter = outIter / 4; if (statIter == 0) statIter = 1;
    iniIter = iniTime = 0;
    maxIter = util::roundToInt(maxTime/(1000.0*delta_t));

    if (useMask) {
        if (!toCluster)
            maskFile = basePath + maskFile;
        plb_ifstream infile(maskFile.c_str());
        if (infile.good()) {
            pcout << "Found mask file!" << std::endl;
            infile >> *mask;
            foundMask = true;
        } else {
            pcout << "Creating mask file..." << std::endl;
        }
    }

    for (plint i = 0; i < objCount; i++) {
        for (plint xyz = 0; xyz < 3; xyz++)
            velocity[i][xyz] *= uLB;
        if (objType[i] == "fluid" || objType[i] == "inflow" || objType[i] == "outflow") {
            bounds[i].x0 = util::roundToInt(xmin[i]*nx/xmax[domainIndex]);
            bounds[i].x1 = util::roundToInt(xmax[i]*nx/xmax[domainIndex]);
            bounds[i].y0 = util::roundToInt(ymin[i]*ny/ymax[domainIndex]);
            bounds[i].y1 = util::roundToInt(ymax[i]*ny/ymax[domainIndex]);
            bounds[i].z0 = util::roundToInt(zmin[i]*nz/zmax[domainIndex]);
            bounds[i].z1 = util::roundToInt(zmax[i]*nz/zmax[domainIndex]);
        }
        if (objType[i] == "fluid") {
            fluidCount++;
            TriangleSet<T> triangleSet(fName[i], DBL);
            triangleSet.translate(Array<T,3>(-dxmin,-dymin,-dzmin));
            triangleSet.scale(lx/dx,ly/dy,lz/dz);
            fluidMesh.push_back(new DEFscaledMesh<T>(triangleSet));
        } else if (objType[i] == "wall") {
            staticSurfaces++;
            TriangleSet<T> triangleSet(fName[i], DBL);
            triangleSet.translate(Array<T,3>(-dxmin,-dymin,-dzmin));
            triangleSet.scale(lx/dx,ly/dy,lz/dz);
            solidMesh.push_back(new DEFscaledMesh<T>(triangleSet));
            std::ostringstream oss;
            oss << "staticobject_" << staticSurfaces-1 << "_";
            std::string fn = oss.str();
            triangleSet.writeBinarySTL(outDir + createFileName(fn, 0, PADDING) + ".stl");
        } else if (objType[i] == "obstacle") {
            movingSurfaces++;
            TriangleSet<T> *triangleSet = new TriangleSet<T>(fName[i], DBL);
            triangleSet->translate(Array<T,3>(-dxmin,-dymin,-dzmin));
            triangleSet->scale(lx/dx,ly/dy,lz/dz);
            plint maxRefinements = MAX_REFINE;
            bool succeeded = triangleSet->refineRecursively(delta_x, maxRefinements);
            if (!succeeded) {
                pcout << std::endl;
                pcout << "WARNING: The target maximum triangle edge length " << delta_x
                      << " for the surface " << i << " was not reached after " << maxRefinements
                      << " refinement iterations." << std::endl;
                pcout << std::endl;
            }
            triangleSet->scale(1.0/delta_x);
            ConnectedTriangleSet<T> connectedTriangleSet(*triangleSet);
            plint numVertices = connectedTriangleSet.getNumVertices();
            startIds.push_back(vertexArray.size());
            for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
                vertexArray.push_back(connectedTriangleSet.getVertex(iVertex));
                T area;
                Array<T,3> unitNormal;
                connectedTriangleSet.computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
                areas.push_back(area);
                flags.push_back(i);
            }
            movingObstacles.push_back(connectedTriangleSet);
/*      } else if (objType[i] == "inflow") {  // --- Removed for causing the "squirt bug" - must disable debug mode!
            if (fabs(velocity[i][0]) > fabs(velocity[i][2])) {
                if (fabs(velocity[i][0]) > fabs(velocity[i][1])) {
                    plint avg = (bounds[i].x0 + bounds[i].x1) / 2;
                    bounds[i].x0 = bounds[i].x1 = avg;
                } else {
                    plint avg = (bounds[i].y0 + bounds[i].y1) / 2;
                    bounds[i].y0 = bounds[i].y1 = avg;
                }
            } else if (fabs(velocity[i][2]) > fabs(velocity[i][1])) {
                plint avg = (bounds[i].z0 + bounds[i].z1) / 2;
                bounds[i].z0 = bounds[i].z1 = avg;
            } else {
                plint avg = (bounds[i].y0 + bounds[i].y1) / 2;
                bounds[i].y0 = bounds[i].y1 = avg;
            } */
        }
    }
    
    if (disableST) {
        surfaceTensionLB = 0.0;
        contactAngle = -0.01;
    }

    // Construct simulation cell

#ifdef PLB_DEBUG
    pcout << "Constructing simulation cell..." << std::endl;
#endif
    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz));
    Dynamics<T,DESCRIPTOR>* dynamics
//      = new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(omega, cSmago);
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, cSmago);
//      = new IncBGKdynamics<T,DESCRIPTOR>(omega);
    FreeSurfaceFields3D<T,DESCRIPTOR>* fields;
    if (movingSurfaces) {
        fields = new FreeSurfaceFields3D<T,DESCRIPTOR>(blockStructure, dynamics->clone(), rhoEmpty,
                surfaceTensionLB, contactAngle, externalForce, ibIter, vertexArray, areas, flags,
                ForwardVelocity(vertexArray), false, 0);
    } else {
        fields = new FreeSurfaceFields3D<T,DESCRIPTOR>(blockStructure, dynamics->clone(), rhoEmpty,
                surfaceTensionLB, contactAngle, externalForce);
    }
    if (useMask && !foundMask) {
        outMask = new int**[nx];
        for (plint xx = 0; xx < nx; xx++) {
            outMask[xx] = new int*[ny];
            for (plint yy = 0; yy < ny; yy++) {
                outMask[xx][yy] = new int[nz];
                for (plint zz = 0; zz < nz; zz++)
                    outMask[xx][yy][zz] = freeSurfaceFlag::wall;
            }
        }
    }
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    boundaryCondition->setVelocityConditionOnBlockBoundaries(fields->lattice, fields->lattice.getBoundingBox(),
        boundary::outflow);
    if (foundMask) {
        setToConstant(fields->flag, *mask, (int)freeSurfaceFlag::wall, fields->flag.getBoundingBox(), (int)freeSurfaceFlag::wall);
        setToConstant(fields->flag, *mask, (int)freeSurfaceFlag::fluid, fields->flag.getBoundingBox(), (int)freeSurfaceFlag::fluid);
        setToConstant(fields->flag, *mask, (int)freeSurfaceFlag::interface, fields->flag.getBoundingBox(), (int)freeSurfaceFlag::interface);
        setToConstant(fields->flag, *mask, (int)freeSurfaceFlag::empty, fields->flag.getBoundingBox(), (int)freeSurfaceFlag::empty);
    } else {
        setToConstant(fields->flag, fields->flag.getBoundingBox(), (int)freeSurfaceFlag::wall);
        setToFunction(fields->flag, fields->flag.getBoundingBox().enlarge(-1), initialFluidFlags);
//      setToFunction(fields->flag, fields->flag.getBoundingBox(), initialFluidFlags);
    }
    if (isInlet) setToConstant(fields->flag, bounds[inletIndex], (int)freeSurfaceFlag::fluid);

    fields->lattice.toggleInternalStatistics(false);
    fields->periodicityToggleAll(false);
    fields->defaultInitialize();

    if (useMask) {
        if (!foundMask) {
            plb_ofstream omf(maskFile.c_str());
            for (plint xx = 0; xx < nx; xx++)
                for (plint yy = 0; yy < ny; yy++)
                    for (plint zz = 0; zz < nz; zz++)
                        omf << fields->flag.get(xx,yy,zz) << " ";
            global::mpi().barrier();
        }

        delete mask;
        if (!foundMask) {
            for (plint xx = 0; xx < nx; xx++) {
                for (plint yy = 0; yy < ny; yy++)
                    delete outMask[xx][yy];
                delete outMask[xx];
            }
            delete outMask;
        }
    }

    // Load the previous final state for continuous simulation (or not)

    if (continuous) {
        pcout << "Checking for saved states" << std::endl;
        // TODO: parse numeric suffix using readlink() from linked filename to iniIter
        char buf[100];
        memset(buf, 0, sizeof(buf));
        if (readlink(LAT_SAVE, buf, sizeof(buf)-1) > 0) {
            pcout << "State found! Loading..." << std::endl;
            loadBinaryBlock(fields->flag, FLAG_SAVE);
            loadBinaryBlock(fields->j, J_SAVE);
            loadBinaryBlock(fields->lattice, LAT_SAVE);
            loadBinaryBlock(fields->mass, MASS_SAVE);
            loadBinaryBlock(fields->normal, NORM_SAVE);
            loadBinaryBlock(fields->rhoBar, RHO_SAVE);
            loadBinaryBlock(fields->volumeFraction, VF_SAVE);
            for (plint i = 8; isdigit(buf[i]); i++)
                iniTime = iniTime*10 + buf[i] - '0';
            maxTime += iniTime;
            iniIter = (int)(iniTime/(1000.0*delta_t));
            maxIter += iniIter;
        } else {
            pcout << "State not found. Starting new simulation..." << std::endl;
        }
    }

    // Initialization for fluid movement and inflow/outflow

#ifdef PLB_DEBUG
    pcout << "Initializing fluids..." << std::endl;
#endif
    for (plint i = 0; i < objCount; i++) {
        if ((objType[i] == "fluid" || objType[i] == "inflow")
                    && !(velocity[i][0] == 0.0 && velocity[i][1] == 0.0 && velocity[i][2] == 0.0)) {
            boundaryCondition->addVelocityBoundary0N(bounds[i], fields->lattice); // Check normal orientation!
            setBoundaryVelocity(fields->lattice, bounds[i], velocity[i]);
        }
    }
    if (isOutlet) defineDynamics(fields->lattice, bounds[outletIndex], new BoundaryCompositeDynamics<T,DESCRIPTOR>(dynamics->clone(), false) );

    // Main iteration loop

    pcout << "Starting simulation..." << std::endl;
    saveIndex = iniIter/outIter + saveOffset;
    isInlet = checkInlet(saveIndex, enabledFrame, enabledSwitch);
    rotNorm = FPS * delta_t;
    velNorm = rotNorm * (lx/dx) / delta_x;
    for (iT = iniIter; iT <= maxIter; iT++) {
        T sum_of_mass_matrix = T();
        T lost_mass = T();
        if (iT % statIter == 0 || iT == maxIter) {
            pcout << "Simulation time: " << iT*1000.0*delta_t << "/" << maxTime << " milliseconds" << std::endl;
            fields->lattice.toggleInternalStatistics(true);
        }
        bool tempInlet = isInlet;
        if (iT % outIter == 0 || iT == maxIter) {
            pcout << "Writing results at time " << iT * delta_t << " sec (frame " << saveIndex << ")" << std::endl;
            writeResults(fields, iT);
            saveIndex++;
            isInlet = checkInlet(saveIndex, enabledFrame, enabledSwitch);
            if (STLmov) {
                for (plint iSurface = 0; iSurface < movingSurfaces; iSurface++) {
                    TriangleSet<T> *triangleSet = movingObstacles[iSurface].toTriangleSet(DBL, &vertexArray,
                            startIds[iSurface]);
                    PLB_ASSERT(triangleSet != 0);
                    triangleSet->scale(delta_x);
                    std::ostringstream oss;
                    oss << "movingobject_" << iSurface << "_";
                    std::string fn = oss.str();
                    triangleSet->writeBinarySTL(outDir + createFileName(fn, (int)(iT*1000.0*delta_t), PADDING) + "ms.stl");
                    delete triangleSet;
                }
            }
        }

        fields->lattice.executeInternalProcessors();
        fields->lattice.incrementTime();

        if (iT % statIter == 0 || iT == maxIter) {
            fields->lattice.evaluateStatistics();
            sum_of_mass_matrix = fields->lattice.getInternalStatistics().getSum(0);
            fields->lattice.toggleInternalStatistics(false);
            if (isnan(sum_of_mass_matrix)) {
                pcout << "Simulation crashed! Exiting..." << std::endl;
                exit(1);
            }
        }
        if (isInlet) {
            if (!tempInlet) {
                setToConstant(fields->flag, bounds[inletIndex], (int)freeSurfaceFlag::fluid);
                applyProcessingFunctional(new PartiallyDefaultInitializeFreeSurface3D<T,DESCRIPTOR>(dynamics->clone(), externalForce, rhoEmpty),
                    bounds[inletIndex].enlarge(1), fields->freeSurfaceArgs);
                OnLatticeBoundaryCondition3D<T,DESCRIPTOR>*
                    newBoundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
                newBoundaryCondition->addVelocityBoundary0N(bounds[inletIndex], fields->lattice);
                setBoundaryVelocity(fields->lattice, bounds[inletIndex], velocity[inletIndex]);
                delete newBoundaryCondition;
            }
            applyProcessingFunctional(new InletConstVolumeFraction3D<T,DESCRIPTOR>(1.002), bounds[inletIndex], fields->freeSurfaceArgs);
        }
        if (isOutlet)
            applyProcessingFunctional(new OutletMaximumVolumeFraction3D<T,DESCRIPTOR>(0.5), bounds[outletIndex], fields->freeSurfaceArgs);

        // Immersed Boundary Method.

        plint nextIt = iT + 1;

        for (plint iSurface = 0; iSurface < movingSurfaces; iSurface++) {
            plint numVertices = (plint) movingObstacles[iSurface].getNumVertices();
            plint startId = startIds[iSurface];
            for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
                plint id = iVertex + startId;
                moveStep = calculateVelocity(saveIndex-1, vertexArray[id], moveX, moveY, moveZ, rotX, rotY, rotZ);
                vertexArray[id] += moveStep;
            }
        }

    }

    pcout << "Simulation finished!" << std::endl;

    // Save the final state as checkpoint for continuous simulation

    if (continuous) {
        pcout << "Saving simulation state" << std::endl;

        sstm << outDir << "flag_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->flag, dat);

        sstm << outDir << "j_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->j, dat);

        sstm << outDir << "lattice_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->lattice, dat);

        sstm << outDir << "mass_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->mass, dat);

        sstm << outDir << "normal_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->normal, dat);

        sstm << outDir << "rhoBar_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->rhoBar, dat);

        sstm << outDir << "volFrac_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        saveBinaryBlock(fields->volumeFraction, dat);

        // Make symbolic links for easier access when loading

        remove(FLAG_SAVE);
        sstm << "flag_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), FLAG_SAVE);

        remove(J_SAVE);
        sstm << "j_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), J_SAVE);

        remove(LAT_SAVE);
        sstm << "lattice_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), LAT_SAVE);

        remove(MASS_SAVE);
        sstm << "mass_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), MASS_SAVE);

        remove(NORM_SAVE);
        sstm << "normal_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), NORM_SAVE);

        remove(RHO_SAVE);
        sstm << "rhoBar_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), RHO_SAVE);

        remove(VF_SAVE);
        sstm << "volFrac_" << maxTime;
        dat = sstm.str();
        sstm.str(""); sstm.clear();
        symlink(dat.c_str(), VF_SAVE);
    }

    delete boundaryCondition;
    delete dynamics;
    delete fields;

    exit(0);

}
