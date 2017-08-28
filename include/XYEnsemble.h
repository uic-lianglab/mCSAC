/*
 * XYEnsemble.h
 * Author: Yun Xu
 * Email: yxu7@uic.edu
 * Date: Sep 20, 2011
 */
#ifndef XYENSEMBLE_H
#define XYENSEMBLE_H

#include "XYMatrix.h"
#include "tree.hh"
#include "XYMath.h"
#include "MersenneTwister.h"
#include "XYSO3Sequence.h"
#include <sstream>
#include "octree.h"
#include "libgen.h"

using namespace std;
using namespace spatialaggregate;

typedef pair<float,int> FloatIntPair;
struct FloatIntPairCompare{
  bool operator () ( const FloatIntPair& left, const FloatIntPair& right)
  { return left.first < right.first; }
};

//template<class T1, class T2 >
//struct sort_pair_first_greater {
//	bool operator()(const pair<T1,T2>&left, const pair<T1,T2>&right) {
//		return left.first > right.first;
//	}
//};


class CXYEnsemble{
public:
	CXYEnsemble();
	CXYEnsemble(
		char* cOutPath,
		float fPersistenLength,
		float fCollisionLength,
		float fPackingDensity,
		float fNucleusSphereDiameter,
		int iNumNodes,
		int iNumSamplePoints,
		int NumberofChains,
		vector <int> ChainLengths,
		char* cStartEndFile,
		char* cContIndFile
		);

	~CXYEnsemble(void);

	// set and get bindingangle begin
	void SetBendingAngleBeg(float fAng);
	float GetBendingAngleBeg(void);

	// set and get bindingangle end
	void SetBendingAngleEnd(float fAng);
	float GetBendingAngleEnd(void);

	// set and get number of binding angles
	void SetNumBendingAngle(int iNum);
	int GetNumBendingAngle(void);

	// set and get number of torsion angles
	void SetNumTorsionAngle(int iNum);
	int GetNumTorsionAngle(void);

	// set and get torsionangle begin
	void SetTorsionAngleBeg(float fAng);
	float GetTorsionAngleBeg(void);

	// set and get torsionangle end
	void SetTorsionAngleEnd(float fAng);
	float GetTorsionAngleEnd(void);

	// set and get persistence length
	void SetPersistenceLength(float fPL);
	float GetPersistenceLength(void);

	// set and get collision length
	void SetCollisionLength(float fCollision);
	float GetCollisionLength(void);

	void SetPackingDensity(float fPackingDensity);
	float GetPackingDensity(void);

	// set and get number of nodes
	void SetNumNodes(int iNumNodes);
	int GetNumNodes(void);

	// set and get number of chains
	void SetNumChains(int NumberofChains);
	int GetNumChains(void);

	// set and get number of nodes
	void SetChainLengths(vector <int> ChainLengths);
	vector <int> GetChainLengths(void);



	void SetSamplesOrg();
	CXYMatrix<float> GetSamplesOrg();

	// use V1->V2 as Z axis, reconstruct new xyz coordinate system
	CXYMatrix<float> NewXYZ(CXYVector<float> kV_1, CXYVector<float> kV_2);
	CXYMatrix<float> NewXYZ(CXYVector<float> kV_0, CXYVector<float> kV_1, CXYVector<float> kV_2);
	// get rotation matrix
	CXYMatrix<float> GetRotMatrix(CXYMatrix<float> &rkMXYZ);

	// get node samples based on current node
	CXYMatrix<float> GetNodeSamples(tree<CXYVector<float> >::iterator  &itNode, int iSegInd);


	// initialize chain (0, 0, 0) and (0, 0, PersistenceLength)
	void InitializeChain(int m_samples);


	vector <vector<tree< CXYVector<float> >*> > GetTree();
	CXYVector<float> RndSetStartPoint(CXYVector<float>  kV_centrpoint);
	CXYVector<float> RndSetEndPoint(int j, CXYVector <float> cpoint);
	CXYVector<float> RndSetCentrPoint(int j);
	CXYVector<float> RndSetCentrPoint_nucleolus(void);
	CXYVector <int> getCentromeres(int NofC);
	float GetSBPSphereDiameter(void);
	bool IsInsideSBPSphere(float* fCoord);
	bool IsEqualityHolds(CXYVector<float> point, int j);


	void SetNucleusSphereDiameter(float fNucleusSphereDiameter);
	float GetNucleusSphereDiameter(void);

	bool IsInsideSphere(CXYVector<float> kV_point);
	bool IsInsideSphere(float* fCoord);
	bool IsCollision(CXYVector<float> kV_point, int j);
  bool IsCollision(CXYMatrix<float> & MiddleEndPoints, int j);

  bool GrowOneChain(int m_samples, int j, int numnodes, int a, int dir, CXYVector<float> points);
	// growth chain
 bool GrowmChain(int k, int f, int m_samples);

  int getMaxofArray(CXYVector<int> array_, int length_);

	void Resampling(int k, int m_samples);
	 double CalculateEffSampleSize(int m);

	void WritemChain(char* fn, int m, int m_samples);
	void WriteDistance(char* fn);
	void SetOutPath(char* cPathName);
	char* GetOutPath(void);
	// void SetSegLengths(char* cStartEndFile);
	void SetSegLengths(char* cStartEndFile,const char* cMethod);
	void SetSegLengths(void);
	vector<float> & GetSegLengths(void);
	float GetSegLength(int ind);

	void SetContIndex(char* cContIndFile);
  void SetContIndex(void);

  void GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints);
  bool IsSatisfyCondition(CXYMatrix<float> & MiddleEndPoints, int j);
  bool IsOnNE(float* fCoord);
  bool IsSatisfyCentr(CXYVector<float> & kV_Point, int j, int dir, int i);
  bool IsCloseSBP(CXYVector<float> & kV_Point, float SBPradius, CXYVector<int> & SBPcenter, int centr, int i);
  bool IsCloseNE(CXYVector<float> & kV_Point, int chainlength, int i);
  CXYVector <float> CalculateProbForGrowth(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int excelpointsize, CXYVector <int> excelpointind, int m_samples, int dir, int ind, CXYVector<float> points);
  float CalculateTargetDistribution(CXYVector<float> kV_Point, int dir, CXYVector<float> points);
  float getMinofArray(CXYVector<float> array_, int length_);
  CXYVector <float>  CalculateBeta(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int GoodPointSize, vector<int> goodpointind, int chainid, int dir, int nodeid, CXYVector <float> cent);
  CXYVector <int> getExcellentIndexes(CXYVector <float> beta, int j, int GoodPointSize, vector <int> GoodPointInd, int excellentpoints);
  int getExcellentPoints(CXYVector <float> beta, int j, int GoodPointSize);
  float getNewAngle(CXYVector <float> kV_Point, int chainid, int ind);
  int getExcellentPoints_nucleolus(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int GoodPointSize, vector<int> goodpointind);
  CXYVector <int> getExcellentIndexes_nucleolus(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int GoodPointSize, vector<int> goodpointind, int excellentpoints);


  void GetGoodPoints( CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints,vector<int>& GoodPointInd, vector<int>& NoCollisionPointInd, int j, int dir, int i);

  void CreateNewPointTree(int k, int m_samples, CXYVector<int> vectorofI, int dir);
  void CreateNewOctree(int k, int m_samples);
  CXYMatrix <float> GetPoints(tree <CXYVector <float> >* ptr, int k);
  void ReconstructLogweight(int m_samples, CXYVector<int> vectorofI);
  CXYVector <int> getResampleIndex(int m_samples);
  CXYVector <float> getResampleprob(int m_samples);
  float getInitialAngle(CXYVector <float> cpoint, int j);

  void WriteWeight(char* fn);
private:
	char* m_cOutPath;
	float m_fBendingAngleBeg; // binding angle begin
	float m_fBendingAngleEnd; // binding angle end
	float m_fTorsionAngleBeg;	// torsion angle begin
	float m_fTorsionAngleEnd;	// torsion angle end
	int   m_iNumBendingAngle; 	// number of binding angle
	int   m_iNumTorsionAngle; 	// number of torsion angle
	float m_fPersistenceLength; // persistence length
	float m_fCollisionLength;   // collision length
	float m_fPackingDensity;   // packing density
	float m_fNucleusSphereDiameter;
	int   m_iNumNodes;
	int   m_iNumSamplePoints;
	int m_fNumofChains;
	vector <int> m_fChainLengths;
	CXYMatrix<float>* m_pMSamplesOrg;
	float endpoints[1][32][3];
    float centromeres[1][32][3];


	vector < vector <tree <CXYVector <float> >* > >  m_pTChain; // tree stored coordinates
    float alpha[32];



	vector<float> m_vfSegLength; // vector of segment length (may different)
	CXYMatrix<float> m_MContInd; // vector of contact index

   vector < OcTree<float,int>* > m_pOctree;
  Eigen::Matrix< float, 4, 1 > m_center;
  float m_minimumVolumeSize;
  float m_dr;
  int m_maxDepth;

  int m_iNumMiddleEndPoints;
  vector <FloatIntPair> m_LogWeight;
  vector <FloatIntPair> ms_logweight;

  vector < double> m_vLogWeight;
  int mNum;
  int counter;
  int indicator;
};

#endif
