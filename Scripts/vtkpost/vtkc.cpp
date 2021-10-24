#include <vtkVersion.h>
#include <vtkVersionMacros.h>
#include <vtkSmartPointer.h>
#include <vtkMPIController.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPMultiBlockDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkThreshold.h>
#include <vector>
#include <stdio.h>
#include <typeinfo> // For std::bad_cast
#include <iostream> // For std::cout, std::err, std::endl etc.
#include <fstream>
using namespace std;

class CVTKInterface {

	public:
		CVTKInterface ();
		~CVTKInterface ();
		void VTK_BeginBlocks(int a);
		void VTK_SetBlockName(int a, char *name);
		void VTK_EndBlocks();
		void VTK_BeginPoints(int a);
		void VTK_WritePoints (int a, double x, double y, double z);
		void VTK_EndPoints ();
		void VTK_BeginMesh (int a);
		void VTK_SetCellType (int a);
		void VTK_BeginCell (int a);
		void VTK_WriteCells (int a);
		void VTK_EndCells (int a, int b);
		void initGrid();
		void initBlocks();
		void VTK_BeginWriter (int a, int wType);
		void VTK_BeginParallelWriter (int a);
		void VTK_BeginBlockWriter (int a);
		void VTK_SetWriteFile (char *name);
		void VTK_WriteMesh ();
		void VTK_Flush () {Writer->Write();}
		void VTK_SetMPI (int a, int b);
                void VTK_SetMulticomm (int kfl_multicomm, int MulticommColor) ;
		void VTK_Reset ();
		void VTK_SetTimeWriter (char *name);
		void VTK_AddTime (double x,int a,char *name);
		void VTK_EndTimeWriter (char *name);
		void VTK_AddCycle (int a);
		void VTK_BeginScalar (char *name);
		void VTK_WriteScalar (int a, double x);
		void VTK_EndScalar ();
		void VTK_BeginVector (char *name);
		void VTK_WriteVector (int a, double x, double y, double z);
		void VTK_EndVector ();
		void VTK_BeginArray4 (char *name);
		void VTK_WriteArray4 (int a, double x, double y, double z, double w);
		void VTK_EndArray4 ();
		void VTK_BeginArray6 (char *name);
		void VTK_WriteArray6 (int a, double x, double y, double z, double w, double v, double t);
		void VTK_EndArray6 ();
		void VTK_BeginArray9 (char *name);
		void VTK_WriteArray9 (int a, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9);
		void VTK_EndArray9 ();
		void VTK_BeginField (char *name);
		void VTK_WriteField (const double *x);
		void VTK_EndField ();
		void VTK_BeginScalarGP (char *name);
		void VTK_WriteScalarGP (int a, double x);
		void VTK_EndScalarGP ();
		void VTK_BeginVectorGP (char *name);
		void VTK_WriteVectorGP (int a, double x, double y, double z);
		void VTK_EndVectorGP ();
		void VTK_BeginArray4GP (char *name);
		void VTK_WriteArray4GP (int a, double x, double y, double z, double w);
		void VTK_EndArray4GP ();
		void VTK_BeginArray6GP (char *name);
		void VTK_WriteArray6GP (int a, double x, double y, double z, double w, double v, double t);
		void VTK_EndArray6GP ();
		void VTK_BeginArray9GP (char *name);
		void VTK_WriteArray9GP (int a, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9);
		void VTK_EndArray9GP ();

		//Getters
		vtkSmartPointer<vtkXMLWriter> GetParallelWriter();

	private:

		//General variables
		vtkMPIController* controller;
                vtkMPICommunicator* communicator;
		vtkSmartPointer<vtkMultiBlockDataSet> blocks;
		vtkSmartPointer<vtkPoints> points;
		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
		vtkSmartPointer<vtkUnsignedCharArray> ghostArray;
		vtkSmartPointer<vtkDoubleArray> scalar;
		vtkSmartPointer<vtkDoubleArray> field;
		vtkSmartPointer<vtkDoubleArray> time;
		vtkSmartPointer<vtkIntArray> cycle;
		vtkSmartPointer<vtkDoubleArray> vec;
		vector<int> GlobalIds;
		vector<string> bname;
		//Writers
		int writerType;
		vtkSmartPointer<vtkXMLWriter> Writer;
		vtkSmartPointer<vtkXMLPMultiBlockDataWriter> blockWriter;
		vtkSmartPointer<vtkXMLPUnstructuredGridWriter> parallelWriter;
		//Filters
		vector<vtkSmartPointer<vtkThreshold> > threshold;
		//Extracted mesh
		vector<vtkSmartPointer<vtkUnstructuredGrid> > blockMesh;

		enum { Ascii, Binary, Appended };
		enum { parallel, block};
		int numPoints;
		int numCells;
		int MPIsize;
		int MPIrank;
		int numId;
		int count;
		int CellNumPoints;
		int CellTypeNumerator;
		int argc;
		int numblocks;
		char ** argv;
		string filename;
		string fname;
		string tname;
		ofstream pvdfile;
};

CVTKInterface::CVTKInterface (){
}
CVTKInterface::~CVTKInterface (){
}

void CVTKInterface::VTK_BeginBlocks(int numB) {

	numblocks=numB;
	blocks -> SetNumberOfBlocks(numblocks);

	threshold.resize(numblocks);
	for(int i=0; i<numblocks;++i) {
		threshold[i] = vtkSmartPointer<vtkThreshold>::New();
	}
	blockMesh.resize(numblocks);
}

void CVTKInterface::VTK_SetBlockName(int iblock,char * name) {

	bname[iblock] = std::string(name);
}

void CVTKInterface::VTK_EndBlocks() {

	for(int i=0; i<numblocks;++i) {
		//(blocks-> GetMetaData(i)).Set(vtkCompositeDataSet.NAME(),bname[i]);
		threshold[i]->SetInputData(unstructuredGrid);
		threshold[i]->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"blocks");
		threshold[i]->ThresholdBetween(i-0.5,i+0.5);
		threshold[i]->Update();
		blockMesh[i]=threshold[i]->GetOutput();
		blocks -> SetBlock(i,blockMesh[i]);
	}

}

void CVTKInterface::VTK_BeginPoints (int numP) {
	points = vtkSmartPointer<vtkPoints>::New();
	numPoints = numP;
	points->SetNumberOfPoints(numPoints);
}
void CVTKInterface::VTK_WritePoints (int id, double x, double y, double z) {
	points->SetPoint(id,x,y,z);
}
void CVTKInterface::VTK_EndPoints () {
	unstructuredGrid->SetPoints(points);
}

void CVTKInterface::VTK_BeginMesh (int numelem) {
	numCells = numelem;
	unstructuredGrid->Allocate(numCells);
        ghostArray->SetNumberOfTuples(numCells);
        ghostArray->SetName(vtkDataSetAttributes::GhostArrayName());
}
void CVTKInterface::VTK_SetCellType (int CellType) {
	CellTypeNumerator = CellType;
}
void CVTKInterface::VTK_BeginCell (int numnodes) {
	numId = 0;
	CellNumPoints = numnodes;
	GlobalIds.resize(CellNumPoints);
}
void CVTKInterface::VTK_WriteCells (int gIds) {
	GlobalIds[numId++] = gIds;
}
void CVTKInterface::VTK_EndCells (int cellId, int ghostLevel) {
	vtkIdType ptIds[CellNumPoints];
	for (int i=0; i<CellNumPoints; i++){
		ptIds[i] = GlobalIds[i];
	}
	unstructuredGrid->InsertNextCell(CellTypeNumerator,CellNumPoints,ptIds);
        
        ghostArray->SetValue(cellId,ghostLevel);
}

void CVTKInterface::VTK_SetMPI (int mpisize, int mpirank) {

	controller = vtkMPIController::New();
	controller->Initialize();
	controller->SetGlobalController(controller);

	MPIsize = mpisize;
	MPIrank = mpirank;
}
void CVTKInterface::initGrid() {

	unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	ghostArray = vtkSmartPointer<vtkUnsignedCharArray>::New();

}
void CVTKInterface::initBlocks() {

	blocks = vtkSmartPointer<vtkMultiBlockDataSet>::New();

}
void CVTKInterface::VTK_BeginParallelWriter (int VTKformat) {

	parallelWriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
	parallelWriter->SetNumberOfPieces(MPIsize);
	parallelWriter->SetGhostLevel(1);

	parallelWriter->SetStartPiece(MPIrank);
	parallelWriter->SetEndPiece(MPIrank);

}

void CVTKInterface::VTK_BeginBlockWriter (int VTKformat) {

	blockWriter = vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
	blockWriter->SetNumberOfPieces(MPIsize);
	blockWriter->SetGhostLevel(1);

}

void CVTKInterface::VTK_BeginWriter (int VTKformat,int wType) {

	writerType = wType;
	initGrid();

	switch(writerType){
		case parallel:
			VTK_BeginParallelWriter(VTKformat);
			break;
		case block:
	                initBlocks();
			VTK_BeginBlockWriter(VTKformat);
			break;
	}
	Writer = GetParallelWriter();

	Writer->GetDefaultFileExtension();

	switch(VTKformat){
		case Ascii:
			Writer->SetDataModeToAscii ();
			break;
		case Binary:
			Writer->SetDataModeToBinary ();
			break;
		case Appended:
			Writer->SetDataModeToAppended ();
			break;
	}
}
vtkSmartPointer<vtkXMLWriter> CVTKInterface::GetParallelWriter() {

	switch(writerType){
		case parallel:
			return this->parallelWriter;
			break;
		case block:
			return this->blockWriter;
			break;
	}

}
void CVTKInterface::VTK_WriteMesh () {
	unstructuredGrid->GetCellData()->AddArray(ghostArray);
#if VTK_MAJOR_VERSION <= 5
	switch(writerType){
		case parallel:
			Writer->SetInput(unstructuredGrid);
			break;
		case block:
			Writer->SetInput(blocks);
			break;
	}
#else
	switch(writerType){
		case parallel:
			Writer->SetInputData(unstructuredGrid);
			break;
		case block:
			Writer->SetInputData(blocks);
			break;
	}
#endif
}
void CVTKInterface::VTK_SetWriteFile (char *filename) {
	fname = std::string(filename);
	fname.append(".pvtu");
	Writer->SetFileName(fname.c_str());
}

void CVTKInterface::VTK_SetMulticomm (int kfl_multicomm, int MulticommColor) {
   communicator = vtkMPICommunicator::New();
   communicator = vtkMPICommunicator::GetWorldCommunicator();
   communicator->SplitInitialize(communicator,MulticommColor,MPIrank);
   //controller->Finalize();
   controller = vtkMPIController::New();
   controller->SetCommunicator(communicator);
   //controller->Initialize();
   controller->SetGlobalController(controller);
}

void CVTKInterface::VTK_Reset () {
	unstructuredGrid->Initialize();
}

void CVTKInterface::VTK_SetTimeWriter(char *filename)
{
	tname = std::string(filename);
	tname.append(".pvd");
	pvdfile.open (tname.c_str());
	pvdfile << "<?xml version=\"1.0\"?>"<<endl;
	pvdfile << "<VTKFile type=\"Collection\" version=\"0.1\">"<<endl;
	pvdfile << "<Collection>"<<endl;
	pvdfile.close();
}

void CVTKInterface::VTK_AddTime(double timevalue,int istep,char *filename)
{
	tname = std::string(filename);
	tname.append(".pvd");
	pvdfile.open (tname.c_str(), std::ios_base::app);
	pvdfile << "<DataSet timestep=\""<<timevalue<<"\" group=\"\" part=\""<<istep<<"\" file=\""<<fname<<"\"/>"<<endl;
	pvdfile.close();

}

void CVTKInterface::VTK_EndTimeWriter(char *filename)
{
	tname = std::string(filename);
	tname.append(".pvd");
	pvdfile.open (tname.c_str(), std::ios_base::app);
	pvdfile << "</Collection>"<<endl;
	pvdfile << "</VTKFile>"<<endl;
	pvdfile.close();
}

void CVTKInterface::VTK_AddCycle(int cyclevalue)
{
	cycle = vtkSmartPointer<vtkIntArray>::New();
	cycle->SetName("CYCLE");
	cycle->SetNumberOfTuples(1);
	cycle->SetTuple1(0, cyclevalue);
	unstructuredGrid->GetFieldData()->AddArray(cycle);
}

void CVTKInterface::VTK_BeginScalar (char *scname) {
	scalar = vtkSmartPointer<vtkDoubleArray>::New();
	scalar->SetNumberOfComponents(1);
	scalar->SetNumberOfTuples(numPoints);
	scalar->SetName(scname);
}
void CVTKInterface::VTK_WriteScalar (int id, double sc) {
	scalar->SetTuple1(id,sc);
}
void CVTKInterface::VTK_EndScalar () {
	unstructuredGrid->GetPointData()->AddArray(scalar);
}
void CVTKInterface::VTK_BeginVector (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(3);
	vec->SetNumberOfTuples(numPoints);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteVector (int id, double x, double y, double z) {
	vec->SetTuple3(id,x,y,z);
}
void CVTKInterface::VTK_EndVector () {
	unstructuredGrid->GetPointData()->AddArray(vec);
}
void CVTKInterface::VTK_BeginField (char *fieldname) {
	field = vtkSmartPointer<vtkDoubleArray>::New();
	field->SetNumberOfComponents(1);
	field->SetName(fieldname);
}
void CVTKInterface::VTK_WriteField (const double *fd) {
	field->InsertNextTuple(fd);
}
void CVTKInterface::VTK_EndField () {
	unstructuredGrid->GetFieldData()->AddArray(field);
}
void CVTKInterface::VTK_BeginArray4 (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(4);
	vec->SetNumberOfTuples(numPoints);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteArray4 (int id, double x1, double x2, double x3, double x4) {
	vec->InsertTuple4(id,x1,x2,x3,x4);
}
void CVTKInterface::VTK_EndArray4 () {
	unstructuredGrid->GetPointData()->AddArray(vec);
}
void CVTKInterface::VTK_BeginArray6 (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(6);
	vec->SetNumberOfTuples(numPoints);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteArray6 (int id, double x1, double x2, double x3, double x4, double x5, double x6) {
	vec->InsertTuple6(id,x1,x2,x3,x4,x5,x6);
}
void CVTKInterface::VTK_EndArray6 () {
	unstructuredGrid->GetPointData()->AddArray(vec);
}
void CVTKInterface::VTK_BeginArray9 (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(9);
	vec->SetNumberOfTuples(numPoints);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteArray9 (int id, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9) {
	vec->InsertTuple9(id,x1,x2,x3,x4,x5,x6,x7,x8,x9);
}
void CVTKInterface::VTK_EndArray9 () {
	unstructuredGrid->GetPointData()->AddArray(vec);
}

void CVTKInterface::VTK_BeginScalarGP (char *scname) {
	scalar = vtkSmartPointer<vtkDoubleArray>::New();
	scalar->SetNumberOfComponents(1);
	numCells = unstructuredGrid->GetNumberOfCells();
	scalar->SetNumberOfTuples(numCells);
	scalar->SetName(scname);
}
void CVTKInterface::VTK_WriteScalarGP (int id, double sc) {
	scalar->InsertTuple1(id,sc);
}
void CVTKInterface::VTK_EndScalarGP () {
	unstructuredGrid->GetCellData()->AddArray(scalar);
}
void CVTKInterface::VTK_BeginVectorGP (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(3);
	numCells = unstructuredGrid->GetNumberOfCells();
	vec->SetNumberOfTuples(numCells);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteVectorGP (int id, double x, double y, double z) {
	vec->InsertTuple3(id,x,y,z);
}
void CVTKInterface::VTK_EndVectorGP () {
	unstructuredGrid->GetCellData()->AddArray(vec);
}
void CVTKInterface::VTK_BeginArray4GP (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(4);
	numCells = unstructuredGrid->GetNumberOfCells();
	vec->SetNumberOfTuples(numCells);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteArray4GP (int id, double x1, double x2, double x3, double x4) {
	vec->InsertTuple4(id,x1,x2,x3,x4);
}
void CVTKInterface::VTK_EndArray4GP () {
	unstructuredGrid->GetCellData()->AddArray(vec);
}
void CVTKInterface::VTK_BeginArray6GP (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(6);
	numCells = unstructuredGrid->GetNumberOfCells();
	vec->SetNumberOfTuples(numCells);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteArray6GP (int id, double x1, double x2, double x3, double x4, double x5, double x6) {
	vec->InsertTuple6(id,x1,x2,x3,x4,x5,x6);
}
void CVTKInterface::VTK_EndArray6GP () {
	unstructuredGrid->GetCellData()->AddArray(vec);
}
void CVTKInterface::VTK_BeginArray9GP (char *vecname) {
	vec = vtkSmartPointer<vtkDoubleArray>::New();
	vec->SetNumberOfComponents(9);
	numCells = unstructuredGrid->GetNumberOfCells();
	vec->SetNumberOfTuples(numCells);
	vec->SetName(vecname);
}
void CVTKInterface::VTK_WriteArray9GP (int id, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9) {
	vec->InsertTuple9(id,x1,x2,x3,x4,x5,x6,x7,x8,x9);
}
void CVTKInterface::VTK_EndArray9GP () {
	unstructuredGrid->GetCellData()->AddArray(vec);
}

extern "C" {
	CVTKInterface *C_VTKInterface_new () {
		return new CVTKInterface();
	}
	CVTKInterface *C_VTKInterface_des (CVTKInterface *This) {
		delete This;
	}
	void C_VTK_BeginBlocks(CVTKInterface *This, int numblocks) {
		This->VTK_BeginBlocks(numblocks);
	}
	void C_VTK_SetBlockName(CVTKInterface *This, int iblock, int nameL, char * name) {
		char auxname[nameL];
		strcpy(auxname,name);
		This->VTK_SetBlockName(iblock,auxname);
	}
	void C_VTK_EndBlocks(CVTKInterface *This) {
		This->VTK_EndBlocks();
	}
	void C_VTK_BeginPoints (CVTKInterface *This, int numPoints) {
		This->VTK_BeginPoints(numPoints);
	}
	void C_VTK_WritePoints (CVTKInterface *This, int id, double x, double y, double z) {
		This->VTK_WritePoints(id,x,y,z);
	}
	void C_VTK_EndPoints (CVTKInterface *This) {
		This->VTK_EndPoints();
	}
	void C_VTK_BeginMesh (CVTKInterface *This, int numelem) {
		This->VTK_BeginMesh(numelem);
	}
	void C_VTK_SetCellType (CVTKInterface *This, int CellType) {
		This->VTK_SetCellType(CellType);
	}
	void C_VTK_BeginCells (CVTKInterface *This, int numnodes) {
		This->VTK_BeginCell(numnodes);
	}
	void C_VTK_WriteCells (CVTKInterface *This, int gIds) {
		This->VTK_WriteCells(gIds);
	}
	void C_VTK_EndCells (CVTKInterface *This, int cellId, int ghostLevel) {
		This->VTK_EndCells(cellId,ghostLevel);
	}
	void C_VTK_BeginWriter(CVTKInterface *This, int format, int wType) {
		This->VTK_BeginWriter(format,wType);
	}
	void C_VTK_SetWriteFile (CVTKInterface *This, int charlen, char * argv) {
		char Filename[charlen];
		strcpy(Filename,argv);
		This->VTK_SetWriteFile(Filename);
	}
	void C_VTK_WriteMesh (CVTKInterface *This) {
		This->VTK_WriteMesh();
	}
	void C_VTK_Flush (CVTKInterface *This) {
		This->VTK_Flush();
	}
	void C_VTK_SetMPI (CVTKInterface *This, int mpisize, int mpirank) {
		This->VTK_SetMPI(mpisize,mpirank);
	}
   void C_VTK_SetMulticomm (CVTKInterface *This, int kfl_multicomm, int MulticommColor) {
		This->VTK_SetMulticomm(kfl_multicomm,MulticommColor);
	}


	void C_VTK_Reset (CVTKInterface *This) {
		This->VTK_Reset();
	}
	void C_VTK_SetTimeWriter (CVTKInterface *This, int charlen, char * argv) {
		char Filename[charlen];
		strcpy(Filename,argv);
		This->VTK_SetTimeWriter(Filename);
	}
	void C_VTK_AddTime (CVTKInterface *This, int charlen, char * argv, double time, int istep) {
		char Filename[charlen];
		strcpy(Filename,argv);
		This->VTK_AddTime(time,istep,Filename);
	}
	void C_VTK_EndTimeWriter (CVTKInterface *This, int charlen, char * argv) {
		char Filename[charlen];
		strcpy(Filename,argv);
		This->VTK_EndTimeWriter(Filename);
	}
	void C_VTK_AddCycle (CVTKInterface *This, int cycle) {
		This->VTK_AddCycle(cycle);
	}
	void C_VTK_BeginScalar (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginScalar(Scalarname);
	}
	void C_VTK_WriteScalar (CVTKInterface *This, int pointid, double sc) {
		This->VTK_WriteScalar(pointid,sc);
	}
	void C_VTK_EndScalar (CVTKInterface *This) {
		This->VTK_EndScalar();
	}
	void C_VTK_BeginScalarGP (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginScalarGP(Scalarname);
	}
	void C_VTK_WriteScalarGP (CVTKInterface *This, int elemid, double sc) {
		This->VTK_WriteScalarGP(elemid,sc);
	}
	void C_VTK_EndScalarGP (CVTKInterface *This) {
		This->VTK_EndScalarGP();
	}
	void C_VTK_BeginVector (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginVector(Scalarname);
	}
	void C_VTK_WriteVector (CVTKInterface *This, int pointid, double x, double y, double z) {
		This->VTK_WriteVector(pointid,x,y,z);
	}
	void C_VTK_EndVector (CVTKInterface *This) {
		This->VTK_EndVector();
	}
	void C_VTK_BeginArray4 (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginArray4(Scalarname);
	}
	void C_VTK_WriteArray4 (CVTKInterface *This, int pointid, double x, double y, double z, double w) {
		This->VTK_WriteArray4(pointid,x,y,z,w);
	}
	void C_VTK_EndArray4 (CVTKInterface *This) {
		This->VTK_EndArray4();
	}
	void C_VTK_BeginArray6 (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginArray6(Scalarname);
	}
	void C_VTK_WriteArray6 (CVTKInterface *This, int pointid, double x, double y, double z, double w, double v, double t) {
		This->VTK_WriteArray6(pointid,x,y,z,w,v,t);
	}
	void C_VTK_EndArray6 (CVTKInterface *This) {
		This->VTK_EndArray6();
	}
	void C_VTK_BeginArray9 (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginArray9(Scalarname);
	}
	void C_VTK_WriteArray9 (CVTKInterface *This, int pointid, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9) {
		This->VTK_WriteArray9(pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9);
	}
	void C_VTK_EndArray9 (CVTKInterface *This) {
		This->VTK_EndArray9();
	}
	void C_VTK_BeginVectorGP (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginVectorGP(Scalarname);
	}
	void C_VTK_WriteVectorGP (CVTKInterface *This, int elemid, double x, double y, double z) {
		This->VTK_WriteVectorGP(elemid,x,y,z);
	}
	void C_VTK_EndVectorGP (CVTKInterface *This) {
		This->VTK_EndVectorGP();
	}
	void C_VTK_BeginArray4GP (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginArray4GP(Scalarname);
	}
	void C_VTK_WriteArray4GP (CVTKInterface *This, int elemid, double x, double y, double z, double w) {
		This->VTK_WriteArray4GP(elemid,x,y,z,w);
	}
	void C_VTK_EndArray4GP (CVTKInterface *This) {
		This->VTK_EndArray4GP();
	}
	void C_VTK_BeginArray6GP (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginArray6GP(Scalarname);
	}
	void C_VTK_WriteArray6GP (CVTKInterface *This, int elemid, double x, double y, double z, double w, double v, double t) {
		This->VTK_WriteArray6GP(elemid,x,y,z,w,v,t);
	}
	void C_VTK_EndArray6GP (CVTKInterface *This) {
		This->VTK_EndArray6GP();
	}
	void C_VTK_BeginArray9GP (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginArray9GP(Scalarname);
	}
	void C_VTK_WriteArray9GP (CVTKInterface *This, int elemid, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9) {
		This->VTK_WriteArray9GP(elemid,x1,x2,x3,x4,x5,x6,x7,x8,x9);
	}
	void C_VTK_EndArray9GP (CVTKInterface *This) {
		This->VTK_EndArray9GP();
	}
	void C_VTK_BeginField (CVTKInterface *This, int charlen, char *argv) {
		char Scalarname[charlen];
		strcpy(Scalarname,argv);
		This->VTK_BeginField(Scalarname);
	}
	void C_VTK_WriteField (CVTKInterface *This, const double *fd) {
		This->VTK_WriteField(fd);
	}
	void C_VTK_EndField (CVTKInterface *This) {
		This->VTK_EndField();
	}
}
