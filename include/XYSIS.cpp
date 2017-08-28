#include "XYSIS.h"

//------------------------------------------------------------------------
CXYSIS::CXYSIS()
{
}

CXYSIS::CXYSIS(
  char* cOutPath,
  float fPersistenLength,
  float fCollisionLength,
  float fCutoff,
  float fPackingDensity,
  float fNucleusSphereDiameter,
  int iNumNodes,
  int iNumSamplePoints,
  char* cStartEndFile,
  int iMmax,
  float fRho_1,
  float fRho_2,
  float fRho_3,
  float fTau_t,
  float fAdjust,
  char* cPvalFile,
  char* cContIndFile)
{
  m_cOutPath = new char[1024];
  m_cDistFileName = new char[1024];
  SetOutPath(cOutPath);
  SetPersistenceLength( fPersistenLength );
  SetPackingDensity(fPackingDensity);
  SetNucleusSphereDiameter(fNucleusSphereDiameter);

  m_iNumSamplePoints = iNumSamplePoints;
  if (strcmp(cStartEndFile, "") == 0){
    // coarse version
    // given number of nodes
    SetNumNodes(iNumNodes); 
    // set each segment length equal to persistence length
    SetSegLengths();
    SetContIndex();  
  } else {
    // fine version
    // read start end position file
    // dynamic setting each segment length determined by mass density
    // SetSegLengths(cStartEndFile,"DIF"); 
    SetSegLengths(cStartEndFile, "AVG");
    SetNumNodes(m_vfSegLength.size()); // set node numbers
    SetContIndex(cContIndFile);  
  }
  SetCollisionLength( fCollisionLength );
  m_pTChain = new tree< CXYVector<float> > ();
  m_pMSamplesOrg = new CXYMatrix<float>;
  SetSamplesOrg();
  
  // For SIS
  m_fInvNumber = iMmax * iNumNodes;
  SetMmax( iMmax );
  SetRho_1( fRho_1 );
  SetRho_2( fRho_2 );
  SetRho_3( fRho_3);
  SetTau_t( fTau_t );
  SetAjust( fAdjust );

//  m_pMDist = new CXYMatrix<float>();
//  SetDistFileName(cDistFile);
//  SetDistMatrix();
  
//  SetSegSegPval(cPvalFile);
  SetSegSegPval_2(cPvalFile);
  m_iNumMiddleEndPoints = max(int(m_fPersistenceLength / m_fCollisionLength) - 1,1);
  
  m_pVErr = new CXYVector<float>(iMmax);
  m_pVEColl = new CXYVector<float>(iMmax);

  // Set octree parameter
  m_center = Eigen::Matrix< float, 4, 1 > (0.0f,0.0f,0.0f,1);
  m_minimumVolumeSize = m_fCollisionLength/4; 
  m_dr = m_fCollisionLength*2;
//  m_dr = fCutoff;
  m_Dr = fCutoff;
//  cout << "m_Dr = " << m_Dr << endl;
  boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =   
  boost::make_shared< OcTreeNodeAllocator< float , int > >();
  OcTree<float,int> octree_tmp(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
  m_maxDepth = (int) ceil(octree_tmp.depthForVolumeSize(m_minimumVolumeSize));
  
}

//------------------------------------------------------------------------
CXYSIS::~CXYSIS(void)
{
  // Record which octree has been deleted so that no need do again
  vector<size_t> vec_deleted; 
  
  vector<size_t>::iterator it;
  for (TreeIterator_VecpOcTree_Map::iterator map_it=m_Treeit_VecpOctree_Map.begin(); map_it !=m_Treeit_VecpOctree_Map.end(); map_it++) {
    for (vector<OcTree<float,int>* >::iterator vec_it=(map_it->second).begin(); vec_it!=(map_it->second).end(); vec_it++) {
      if (find(vec_deleted.begin(), vec_deleted.end(), (size_t)(*vec_it)) == vec_deleted.end()) {
        vec_deleted.push_back((size_t) *vec_it);
//        cout << (size_t) *vec_it << endl;
        delete *vec_it;
      }
    }
  }
  
  delete m_pVErr;
  delete m_pVEColl;
//  delete m_pMDist;
  delete m_cDistFileName;
  delete m_cOutPath;
  delete m_pMSamplesOrg;
  delete m_pTChain;
  
}

//------------------------------------------------------------------------
void CXYSIS::SetPersistenceLength(float fPL)
{
  m_fPersistenceLength = fPL;
}
//------------------------------------------------------------------------
float CXYSIS::GetPersistenceLength(void)
{
  return m_fPersistenceLength;
}
//------------------------------------------------------------------------
void CXYSIS::SetCollisionLength(float fCollision)
{
  vector<float>& vfSegment = GetSegLengths();
  float min_diameter = *(min_element(vfSegment.begin(), vfSegment.end()));
  m_fCollisionLength = min(fCollision,min_diameter);
}
//------------------------------------------------------------------------
float CXYSIS::GetCollisionLength(void)
{
  return m_fCollisionLength;
}
//------------------------------------------------------------------------
void CXYSIS::SetPackingDensity(float fPackingDensity)
{
  m_fPackingDensity = fPackingDensity;
}
//------------------------------------------------------------------------
float CXYSIS::GetPackingDensity(void)
{
  return m_fPackingDensity;
}
//------------------------------------------------------------------------
void CXYSIS::SetNucleusSphereDiameter(float fNucleusSphereDiameter)
{
  m_fNucleusSphereDiameter = fNucleusSphereDiameter;
}
//------------------------------------------------------------------------
float CXYSIS::GetNucleusSphereDiameter(void)
{
  return m_fNucleusSphereDiameter;
}
//------------------------------------------------------------------------
void CXYSIS::SetNumNodes(int iNumNodes){
  m_iNumNodes = iNumNodes;
}
//------------------------------------------------------------------------
int CXYSIS::GetNumNodes(void){
  return m_iNumNodes;
}

//------------------------------------------------------------------------
// Generate sphere sample points with radius persistence length.
void CXYSIS::SetSamplesOrg(void)
{
  CXYSO3Sequence sO3sequence(m_iNumSamplePoints);
  sO3sequence.SetSO3Sequence();
  (*m_pMSamplesOrg) = sO3sequence.GetSO3Sequence();
}
//------------------------------------------------------------------------
CXYMatrix<float> CXYSIS::GetSamplesOrg()
{
  return (*m_pMSamplesOrg);
}

//------------------------------------------------------------------------
//CXYMatrix<float> CXYSIS::NewXYZ(CXYVector<float> kV_0, CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//  CXYMatrix<float> kMXYZ(3,3);
//  CXYVector<float> kV_Z = kV_2 - kV_1;
//  CXYVector<float> kV_01 = kV_1 - kV_0;
//  kV_01.Normalize();
//  if (fabs(kV_Z.Dot(kV_01)) < CXYMath<float>::ZERO_TOLERANCE)
//  { // two line parallel
//    kMXYZ = NewXYZ(kV_1,kV_2);
//  } else {
//    float a = kV_Z[0],  b = kV_Z[1],  c = kV_Z[2];
//    float x = kV_01[0], y = kV_01[1], z = kV_01[2];
//    
//    float fX[3] = { y*c - z*b, z*a - x*c, x*b - y*a };
//    CXYVector<float> kV_X(3,fX);
//    kV_X.Normalize();
//    
//    x = kV_Z[0], y = kV_Z[1], z = kV_Z[2];
//    a = kV_X[0], b = kV_X[1], c = kV_X[2];
//    float fY[3] = { y*c - z*b, z*a - x*c, x*b - y*a };
//    CXYVector<float> kV_Y(3,fY);
//    kV_Y.Normalize();
//    kMXYZ.SetRow(0, kV_X);
//    kMXYZ.SetRow(1, kV_Y);
//    kMXYZ.SetRow(2, kV_Z);
//  }
//  
//  return kMXYZ;
//}
////------------------------------------------------------------------------
//CXYMatrix<float> CXYSIS::NewXYZ(CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//  CXYMatrix<float> kMXYZ(3,3);
//  CXYVector<float> kV_Z = kV_2 - kV_1;
//  kV_Z.Normalize();
//  float a = kV_Z[0], b= kV_Z[1], c = kV_Z[2];
//  
//  float x, y, z;
//  if (c > CXYMath<float>::ZERO_TOLERANCE ) 
//  {
//    x = 1; y = 0; z = -(a*x+b*y)/c;
//  } else if (b > CXYMath<float>::ZERO_TOLERANCE ) 
//  {
//    x = 0; z = 1; y = -(a*x+c*z)/b;
//  } else {
//    y = 1; z = 0; x = -(b*y+c*z)/a;
//  }
//  float fX[3] = {x,y,z};
//  CXYVector<float> kV_X(3,fX);
//  float fY[3] = {b*z-c*y, c*x-a*z, a*y-b*x};
//  CXYVector<float> kV_Y(3,fY);
//  
//  kMXYZ.SetRow(0, kV_X);
//  kMXYZ.SetRow(1, kV_Y);
//  kMXYZ.SetRow(2, kV_Z);
//  
//  
//  return kMXYZ;
//}

//------------------------------------------------------------------------
//CXYMatrix<float> CXYSIS::GetRotMatrix(CXYMatrix<float> &rkMXYZ)
//{
//  CXYVector<float> kV_X = rkMXYZ.GetRow(0);
//  CXYVector<float> kV_Y = rkMXYZ.GetRow(1);
//  CXYVector<float> kV_Z = rkMXYZ.GetRow(2);
//  
//  CXYMatrix<float> kM_A(3,3);
//  float alpha, beta, gamma;
//  float Z3diff = fabs(1-kV_Z[2]*kV_Z[2]);
//  if ( Z3diff < CXYMath<float>::ZERO_TOLERANCE) 
//  {
//    alpha = acos( min(float(1.0), max(float(-1.0), kV_X[0])));
//    beta  = 0;
//    gamma = 0;
//  } else {    // http://www.macosxguru.net/article.php?story=20040210124637626
//    float Denorm = sqrt(Z3diff);
//    alpha = acos(min(float(1.0), max(float(-1.0), -kV_Z[1]/Denorm)));
//    beta  = acos(min(float(1.0), max(float(-1.0), -kV_Z[2])));
//    gamma = acos(min(float(1.0), max(float(-1.0), -kV_Y[2]/Denorm)));
//  }
//  
//  kM_A[0][0] = cos(gamma) * cos(alpha) - cos(beta) * sin(alpha) * sin(gamma);
//  kM_A[0][1] = cos(gamma) * sin(alpha) + cos(beta) * cos(alpha) * sin(gamma);
//  kM_A[0][2] = sin(gamma) * sin(beta);
//  kM_A[1][0] = -sin(gamma) * cos(alpha) - cos(beta) * sin(alpha) * cos(gamma);
//  kM_A[1][1] = -sin(gamma) * sin(alpha) + cos(beta) * cos(alpha) * cos(gamma);
//  kM_A[1][2] = cos(gamma) * sin(beta);
//  kM_A[2][0] = sin(beta) * sin(alpha);
//  kM_A[2][1] = -sin(beta) * cos(alpha);
//  kM_A[2][2] = cos(beta);
//  
//  kM_A.GetInverse(kM_A);
//  
//  
//  return kM_A;
//}

//------------------------------------------------------------------------
void CXYSIS::InitializeChain()
{
  tree< CXYVector<float> >* ptr = GetTree();
  // erase all nodes if not empty
  if (!ptr->empty()) {
    ptr->clear();
  }

  
  tree< CXYVector<float> >::iterator top, root, pos;
  
  top = ptr->begin();
  
  // first node position
  CXYVector<float> kV_startpoint = RndSetStartPoint();
  root = ptr->insert(top, kV_startpoint);
  
  CXYMatrix<float> kMSamplesOrg = GetSamplesOrg();
  kMSamplesOrg *= GetSegLength(0); // index is 0, first one segment
  int iSampleSize = kMSamplesOrg.GetRows();
  int iRnd;
  CXYVector<float> kV_Point(3);
  // second node position
  while (1) {
    MTRand mtrand;
    iRnd = mtrand.randInt(iSampleSize-1);  // integer in [0,iSampleSize-1]
    kV_Point = kMSamplesOrg.GetRow(iRnd) + kV_startpoint; // second point
    if (IsInsideSphere(kV_Point) && ! IsCollision(root,kV_Point)) {
      pos = ptr->append_child(root, kV_Point);
      break;
    }
  }
}
//------------------------------------------------------------------------
tree< CXYVector<float> >* CXYSIS::GetTree()
{
  return m_pTChain;
}
//------------------------------------------------------------------------
CXYVector<float> CXYSIS::RndSetStartPoint(void)
{
  float fdiameter = GetNucleusSphereDiameter();
  float fradius = fdiameter/2;
  MTRand mtrand;
  
  // Check if the random point located in the neuclus sphere.
  // When one point satisfy the condition, we accept it.
  float fCoord[3];
  while (1) {
    fCoord[0] = mtrand.randExc(fradius);
    fCoord[1] = mtrand.randExc(fradius);
    fCoord[2] = mtrand.randExc(fradius);
    if (IsInsideSphere(fCoord)){
      break;
    }
  }
//  cout << x << ";" << y << ";" << z<<";" <<endl;
//  cout << x*x+y*y+z*z << endl;
//  cout << fradius2 << endl;
  CXYVector<float> kV (3,fCoord);
  return kV;
}
//------------------------------------------------------------------------
bool CXYSIS::IsInsideSphere(CXYVector<float> kV_point)
{
  float fradius = GetNucleusSphereDiameter()/2;
  float fradius2 = fradius*fradius;
  bool flag;
  if (kV_point.SquaredLength() < fradius2) {
    flag = true;
  }else {
    flag = false;
  }
  return flag;
  
}
//------------------------------------------------------------------------
bool CXYSIS::IsInsideSphere(float* fCoord)
{
  float fradius = GetNucleusSphereDiameter()/2;
  float fradius2 = fradius*fradius;
  bool flag;
  if (fCoord[0]*fCoord[0] + fCoord[1]*fCoord[1] + fCoord[2]*fCoord[2]
      < fradius2) {
    flag = true;
  }else {
    flag = false;
  }
  return flag;
}

//------------------------------------------------------------------------
bool CXYSIS::IsCollision(tree<CXYVector<float> >::iterator  &ritNode, CXYVector<float> kV_point)
{
  float fCollisionLength = GetCollisionLength();
  // to speed up the caculation, we calculate power(2)
  float fCollisionLength_Squar = fCollisionLength*fCollisionLength; 
  // enumerate previous nodes
  tree< CXYVector<float> >* ptr = GetTree();
  tree< CXYVector<float> >::iterator node, lastsecond_node;
  
  tree<CXYVector<float> >::iterator pos = ritNode;
  pos = ptr->parent(pos);
  while (ptr->is_valid(pos)) 
  {
    if( (kV_point-(*pos)).SquaredLength() < fCollisionLength_Squar ){
      return true;
    }
    pos = ptr->parent(pos);
  }
  return false;

}

//------------------------------------------------------------------------
bool CXYSIS::GrowthOneChain()
{
  bool Flag = true;
  InitializeChain();
  MTRand mtrand;
  
  int iNumNodes = GetNumNodes() - 2;
  tree< CXYVector<float> >* ptr = GetTree();
  tree< CXYVector<float> >::pre_order_iterator node;


  for (int i =0 ; i < iNumNodes; i++) {
    
    node = ptr->end(); 
    node --; // go to last node

    // find parent node
//    par = ptr->parent(node);
    
    // find grandparent node
//    gra = ptr->parent(par);

//    CXYMatrix<float> kMXYZ(3);
//    if (ptr->is_valid(gra)) {
//      kMXYZ = NewXYZ((*gra), (*par), (*node));
//    }else {
//      kMXYZ = NewXYZ((*par), (*node));
//    }
//    
//    CXYMatrix<float> kMRotate = GetRotMatrix(kMXYZ);
//
//    CXYMatrix<float> kMSamplesOrg = GetSamplesOrg();
//    CXYMatrix<float> kMSamplesPoints = kMSamplesOrg.TimesTranspose(kMRotate);
    CXYMatrix<float> kMSamplesPoints = GetSamplesOrg();
    kMSamplesPoints *= GetSegLength(i+1); // because initialize chain already done
    
    int iSampleSize = kMSamplesPoints.GetRows();
    int iRnd;

    CXYVector<float> kV_Point(3);
    int n_trial = 1000;
    while (n_trial > 0) {
      iRnd = mtrand.randInt(iSampleSize-1);  // integer in [0,iSampleSize-1]
      kV_Point = kMSamplesPoints.GetRow(iRnd) + (*node); // next point position
      // check if inside neuclus sphere and there is no collision
      if (IsInsideSphere(kV_Point) && (! IsCollision(node,kV_Point))) {
        ptr->append_child(node, kV_Point);
        break;
      }
      n_trial --;
    }
    // do not continue our growing chain
    if (n_trial == 0){
      Flag = false;
      break;
    }
  }
  return Flag;
}
//------------------------------------------------------------------------
void CXYSIS::WriteChain(char* fn)
{
//http://stdcxx.apache.org/doc/stdlibug/34-2.html 
  
  tree<CXYVector<float> > * ptr = GetTree();
  tree< CXYVector<float> >::iterator pos = ptr->begin();
  char buffer[1024];
  
  ostream* fp;
  if (strcmp(fn,"")) {
    sprintf(buffer, "%s/%s", GetOutPath(),fn);
    fp = new std::ofstream(buffer);
  }else {
    fp = &std::cout;
  }

  while (ptr->is_valid(pos)) {
    *fp << (*pos)[0] <<"\t" <<(*pos)[1] <<"\t"<<(*pos)[2]<<endl;
    ++pos;
  }

  if (fp != &std::cout)
    delete fp;  
}
//------------------------------------------------------------------------
void CXYSIS::WriteDistance(char* fn)
{
  //http://stdcxx.apache.org/doc/stdlibug/34-2.html 
  
  tree<CXYVector<float> > * ptr = GetTree();
  tree< CXYVector<float> >::iterator root = ptr->begin();
  tree< CXYVector<float> >::fixed_depth_iterator pos1, pos2;
  float deltaX, deltaY, deltaZ, length;
  char buffer[1024];
  
  ostream* fp;
  if (strcmp(fn,"")) {
    sprintf(buffer, "%s/%s", GetOutPath(),fn);
    fp = new std::ofstream(buffer);
  }else {
    fp = &std::cout;
  }

  int iNumNodes = GetNumNodes();
  for (int i = 0; i<iNumNodes-1; i++) {
    pos1 = ptr->begin_fixed(root, i);
    for (int j = i+1; j<iNumNodes; j++) {
      pos2 = ptr->begin_fixed(root, j);
      
      deltaX = (*pos2)[0]-(*pos1)[0];
      deltaY = (*pos2)[1]-(*pos1)[1];
      deltaZ = (*pos2)[2]-(*pos1)[2];
      length = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
      
      *fp << i+1 << "\t" << j+1 << "\t"<< length << endl;
    }
  }
  
  if (fp != &std::cout)
    delete fp;
}
//------------------------------------------------------------------------
void CXYSIS::SetOutPath(char* cPathName)
{
  CXYFile::MakeDirectory(cPathName,0755);
  strcpy(m_cOutPath,cPathName);
}

//------------------------------------------------------------------------
char* CXYSIS::GetOutPath(void)
{
  return m_cOutPath;
}
//------------------------------------------------------------------------
//void CXYSIS::SetSegLengths(char* cStartEndFile){
//  CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);
//  for (int i= 0; i<kMStartEnd.GetRows(); i++) {
//    int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
//    m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
//  }
//
//}
//------------------------------------------------------------------------
void CXYSIS::SetSegLengths(char* cStartEndFile,const char* cMethod){
  CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);
  
  if (strcmp(cMethod, "AVG") == 0) {
    // average different node length
    float len = 0L;
    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      len = len + ( kMStartEnd[i][1] - kMStartEnd[i][0] +1 ) * GetPackingDensity();
    }
    float avglen = len/(float)kMStartEnd.GetRows();
    // set fixed length
    SetPersistenceLength( avglen );

    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      m_vfSegLength.push_back(GetPersistenceLength());
    }
    
  } else if (strcmp(cMethod, "DIF") == 0) {
    // different node length
    for (int i= 0; i<kMStartEnd.GetRows(); i++) {
      int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
      m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
    }
  } else {
    cerr  << "Please type Method" << endl;
    exit(-1);
  }

  
}

//------------------------------------------------------------------------
vector<float> & CXYSIS::GetSegLengths(void){
  return m_vfSegLength;
}
//------------------------------------------------------------------------
float CXYSIS::GetSegLength(int ind){
  return m_vfSegLength[ind];
}
//------------------------------------------------------------------------
void CXYSIS::SetSegLengths(void){
  for (int i= 0; i<GetNumNodes(); i++) {
    m_vfSegLength.push_back(GetPersistenceLength());
  }
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
void CXYSIS::SetMmax(int iMmax)
{
  m_iMmax = iMmax;
}
//------------------------------------------------------------------------
int CXYSIS::GetMmax(void)
{
  return m_iMmax;
}
//------------------------------------------------------------------------
void CXYSIS::SetRho_1(float fRho_1)
{
  m_fRho_1 = fRho_1;
}
//------------------------------------------------------------------------
float CXYSIS::GetRho_1(void)
{
  return  m_fRho_1;
}
//------------------------------------------------------------------------
void CXYSIS::SetRho_2(float fRho_2)
{
  m_fRho_2 = fRho_2;
}
//------------------------------------------------------------------------
float CXYSIS::GetRho_2(void)
{
  return m_fRho_2;
}
//------------------------------------------------------------------------
void CXYSIS::SetRho_3(float fRho_3)
{
  m_fRho_3 = fRho_3;
}
//------------------------------------------------------------------------
float CXYSIS::GetRho_3(void)
{
  return m_fRho_3;
}
//------------------------------------------------------------------------
void CXYSIS::SetTau_t(float fTau_t)
{
  m_fTau_t = fTau_t;
}
//------------------------------------------------------------------------
float CXYSIS::GetTau_t(void)
{
  return m_fTau_t;
}
//------------------------------------------------------------------------
void CXYSIS::SetAjust(float fAjust)
{
  m_fAdjust = fAjust;
}
//------------------------------------------------------------------------
float CXYSIS::GetAjust(void)
{
  return  m_fAdjust;
}





//------------------------------------------------------------------------
void CXYSIS::SIS_Algorithm(void)
{
  
  MTRand mtrand;
  
  // Sequential important sampling algorithm
  SetSamplesOrg();
  InitializeChain_SISRoot();
  tree< CXYVector<float> >* tr = GetTree();
  tree< CXYVector<float> >::iterator pos;

//  CXYMatrix<float>& kMDist = GetDistMatrix();
//  int iNumSeg = kMDist.GetRows();

  //----------------------------------------------------------------
  // Vist each segment

  int iNumSeg = GetNumNodes();
  for (int iSeg=1; iSeg<iNumSeg; iSeg++) {
    // find upper level node
//    int iLevel_upp = iSeg - 1;
    int im_t = 0;
//    int iCount =  GetCountContactNode(iSeg-1);
    int iCount = GetCountContactNodesFromMap(iSeg);
    
    // get all upper nodes
//    vector< tree<CXYVector<float> >::iterator > vitNodes_upp = GetNodeSet(iLevel_upp);
    
//    vector< CXYMatrix<float> > vChains;
    vector< tree< CXYVector<float> >::iterator > vitPrvNodes;
    
    CXYMatrix<float> kMChain;
    
//    cout << "SegInd: " << iSeg  << ", iCount: " << iCount <<  endl;

    //----------------------------------------------------------------
    // freely grow because no connections under distance range
    if (iCount == 0) {
//      cout << iSeg << endl;
      float fNoConMaxDist = GetNoConMaxDist(iSeg);
      vector<tree<CXYVector<float> >::iterator > vtmp;
      tree<CXYVector<float> >::iterator itNode_cur;
      while (!m_qPrvNodes.empty()) {
        itNode_cur = m_qPrvNodes.front();
        m_qPrvNodes.pop_front(); // pop out
        int iNumTrial = max(m_iNumSamplePoints, 100);

        CXYVector<float> kV_PrvConPoint = GetPrvConPosition(itNode_cur,iSeg);
        
        CXYVector<float> kV_NextVect;
        CXYVector<float> kV_Point(3);
        while (iNumTrial > 0) {
          int iRnd = mtrand.randInt(m_iNumSamplePoints-1);
          kV_NextVect = GetNodeNextPosition(itNode_cur,iSeg,iRnd);
          float Point[] = {kV_NextVect[0],kV_NextVect[1],kV_NextVect[2]};
          kV_Point = CXYVector<float>(3,Point);
          
//          if ((kV_Point - kV_PrvConPoint).Length() < fNoConMaxDist && ! IsCollision(itNode_cur, kV_Point)) {
          if ((kV_Point - kV_PrvConPoint).Length() < fNoConMaxDist) {
//            cout << iSeg << " "<< (kV_Point - kV_PrvConPoint).Length() << " " << fNoConMaxDist << endl;
            break;
          }

          iNumTrial --;
        }
//        if (iNumTrial == 0) {
//          cout << "Give up trial" << endl;
//        }
        
        pos = tr->append_child(itNode_cur,kV_NextVect);
        vtmp.push_back(pos);
      }
      // save node to queue
      for (vector<tree<CXYVector<float> >::iterator>::iterator vtmpit = vtmp.begin(); 
           vtmpit < vtmp.end(); vtmpit++) {
        m_qPrvNodes.push_back(*vtmpit);
      }
      continue; // go to next segment, continue to next segment
      
    }
    
    //----------------------------------------------------------------
    // Nodes to save last points as node
    tree<CXYVector<float> > tree_LV;
    tree<CXYVector<float> >::iterator top_LV,root_LV, pos_LV;
    CXYVector<float> vRoot(0); // Nothing, just represent root
    tree_LV.set_head(vRoot);
    root_LV = tree_LV.begin();
    //----------------------------------------------------------------
    // vector to save chains
    vector< vector<tree<CXYVector<float> >::iterator >* > vvChains;
    //----------------------------------------------------------------

    tree<CXYVector<float> >::iterator itNode_cur;

    while (!m_qPrvNodes.empty()) {
      itNode_cur = m_qPrvNodes.front();
      
      m_qPrvNodes.pop_front();
      
      CXYMatrix<float> kMSamples = GetNodeSamples(itNode_cur,iSeg);
      // Conformations for each node
      for (int iSample=0; iSample<kMSamples.GetRows(); iSample++) {
        CXYVector<float> kVPoint = kMSamples.GetRow(iSample);
        // store previous vector, weight we do not care, so I set to 0
        float ftmp[] = {kVPoint[0],kVPoint[1],kVPoint[2], 0 ,(*itNode_cur)[4],(*itNode_cur)[5],(*itNode_cur)[6]};
        CXYVector<float> kV_LV(7,ftmp);
        
        pos_LV = tree_LV.append_child(root_LV, kV_LV);
        vector<tree<CXYVector<float> >::iterator > *pvOneChain = new vector<tree<CXYVector<float> >::iterator >;

        GrowOneChain_ByVector(itNode_cur,pos_LV,iSeg,pvOneChain);
        vvChains.push_back(pvOneChain);
        vitPrvNodes.push_back(itNode_cur);

      }
    }
    
    
//    vector<tree<CXYVector<float> >::iterator > ppp = *(vvChains[0]);
//    tree<CXYVector<float> >::iterator iit=ppp[0];
//    cout << (*iit)[0] << endl;
    
    //----------------------------------------------------------------
    // Sequence Importance Sampling
    int iL_t = vvChains.size();
    if (iL_t <= GetMmax()){
      
      im_t = iL_t;
      
      for (int iSample=0; iSample<im_t; iSample++) {
        // append child node;
        pos = tr->append_child(vitPrvNodes[iSample], *( (*vvChains[iSample])[0]));
        m_qPrvNodes.push_back(pos);
      }
      
    } 
    else {
      im_t = GetMmax();
      
      // if the node no connection with others, randomly select iMmax node to cont.
      cout << iSeg << " / " << iNumSeg << endl;
//      cout << m_fCollisionLength << endl;
      // calculate beta_t

      vector<float> vfBeta_t;
      CalBeta_t_ByVector(vvChains, iSeg, vfBeta_t);
      
      // binary search const c
      float fConst_C = BinSearchConstC(vfBeta_t);
      
      // get cumsum of min(C*Beta_t,1), order is ascend
      vector<float> vfCumSum, vfMinCBeta1;
      float fCumSum = 0, fMinCBeta1;
      
      for (unsigned int iSample = 0; iSample < vvChains.size(); iSample ++) {
        fMinCBeta1 = min(fConst_C * (vfBeta_t[iSample]),1.0f);
        fCumSum += fMinCBeta1;
        
        vfMinCBeta1.push_back(fMinCBeta1);
        vfCumSum.push_back(fCumSum);
      }

      // We can get a real number in the range 0 to 1, excluding 1
      vector<int> vfSampleSelInd;
      float fRand = 0;
      for (int iSameleSel = 1; iSameleSel <= im_t; iSameleSel++) {
        fRand = iSameleSel - mtrand.randExc();
        
        unsigned int iCumSumInd = 0;
        while (fRand >  vfCumSum[iCumSumInd]) {
          if (iCumSumInd == vfCumSum.size()-1) {
            break;
          }
          ++iCumSumInd;
        }
        vfSampleSelInd.push_back(iCumSumInd);
      }
      
      // append selected coordinates on certain node
      CXYVector<float> kVector(7);
      for (int iSample = 0; iSample < im_t; iSample++) {
        int iChainInd = vfSampleSelInd[iSample];
        kVector = *((*vvChains[iChainInd])[0]);
        // recalculate weight
        kVector[3] = kVector[3] - log(vfMinCBeta1[iChainInd]);
        pos = tr->append_child(vitPrvNodes[iChainInd],kVector);
        m_qPrvNodes.push_back(pos);
        
      }
    } // end if .. else ..
    
    //----------------------------------------------------------------
    // delete the memory we new
    for (unsigned int iNumChainsSize=0; iNumChainsSize <vvChains.size(); iNumChainsSize++) {
      delete vvChains[iNumChainsSize];
    }
    //----------------------------------------------------------------
    
  }

  // output result
  //----------------------------------------------------------------
  GetTargetDistribution();
  //  WritePDBArr();
  WritePtsArr();
}


//------------------------------------------------------------------------
void CXYSIS::SetDistFileName(char* cFileName)
{
  strcpy(m_cDistFileName, cFileName);
}
//------------------------------------------------------------------------
char* CXYSIS::GetDistFileName(void)
{
  return m_cDistFileName;
}
//------------------------------------------------------------------------
void CXYSIS::SetDistMatrix(void)
{
  char* FN = GetDistFileName();
  CXYFile kFile;
  vector<CXYTriple<int,int,float> > vtriple_raw = kFile.ReadSparseMatrixToTriple(FN);
  vector<CXYTriple<int,int,float> >::iterator it_trip;
  //  m_listConNodeInd.clear();
  vector<int> vConNodeIndDup;
  
  for (it_trip = vtriple_raw.begin(); it_trip != vtriple_raw.end(); it_trip++) {
    // set C++ index
    (*it_trip).first = (*it_trip).first - 1;
    (*it_trip).second = (*it_trip).second - 1;
    //    m_listConNodeInd.push_back( (*it_trip).first );
    vConNodeIndDup.push_back((*it_trip).first);
  }
  sort(vConNodeIndDup.begin(), vConNodeIndDup.end());
  vector<int>::iterator new_end = unique(vConNodeIndDup.begin(), vConNodeIndDup.end());

  // get unique sorted connection node index C style
  m_vConNodeInd.clear();
  m_vConNodeInd.assign(vConNodeIndDup.begin(),new_end);
  
//  vector<int>::iterator vit;
//  for (vit = m_vConNodeInd.begin(); vit!=m_vConNodeInd.end(); vit++) {
//    cout << *vit << endl;
//  }
  
//  lowind highind
//  example 
//  0 94
//  0 186
//  0 252
//  94 186  

  for (unsigned int i =0 ; i < m_vConNodeInd.size()-1; i++) {
    for (unsigned int j = i+1; j < m_vConNodeInd.size(); j++) {
//      cout << m_vConNodeInd[i] << " " << m_vConNodeInd[j] << endl;
      for (it_trip = vtriple_raw.begin(); it_trip != vtriple_raw.end(); it_trip++) {
        if ( ((*it_trip).first == m_vConNodeInd[i] && (*it_trip).second == m_vConNodeInd[j] ) ||
             ((*it_trip).first == m_vConNodeInd[j] && (*it_trip).second == m_vConNodeInd[i] ) ) {
          m_vtriple.push_back(CXYTriple<int,int,float>(m_vConNodeInd[i],m_vConNodeInd[j],(*it_trip).third));
//          cout << m_vConNodeInd[i] << " " << m_vConNodeInd[j] << " " << (*it_trip).third << endl;
          break;
        }
      }
    }
  }
  
}
//------------------------------------------------------------------------
CXYMatrix<float>& CXYSIS::GetDistMatrix(void)
{
  return *m_pMDist;
}
//------------------------------------------------------------------------
vector<tree<CXYVector<float> >::iterator > CXYSIS::GetNodeSet(int iLevel)
{
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator root = ptr->begin();
  tree<CXYVector<float> >::fixed_depth_iterator pos = ptr->begin_fixed(root,iLevel);
  
  vector< tree<CXYVector<float> >::iterator > vit;
  while (ptr->is_valid(pos)) 
  {
    //    cout << (*(pos))[3] << endl;
    vit.push_back(pos); 
    ++pos;
  }
  return vit;
}
//------------------------------------------------------------------------
CXYMatrix<float> CXYSIS::GrowthChain(tree< CXYVector<float> >::iterator &ritNode, 
                                     CXYVector<float> &rkVPoint)
{
  cout << sizeof(ritNode) << endl;
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = ritNode;
  
  int iSize = ptr->depth(pos) + 2;
  CXYMatrix<float> kM(iSize,6);
  
  // from bottom to top, fill chain 
  int iRowInd = iSize - 1;
  kM.SetRow(iRowInd, rkVPoint);
  
  while (ptr->is_valid(pos)) 
  {
    iRowInd --;
    kM.SetRow(iRowInd, (*pos));
    pos = ptr->parent(pos);
  }
  
  // according previous node h1 and h2, get new h1 and h2
  kM[iSize-1][4] = h1Function(kM);
  kM[iSize-1][5] = h2Function(kM);
  
  //  char buffer[] = "test.mat";
//  CXYFile::WriteMatrix("test.mat","w",kM);
  
  return kM;
}
//------------------------------------------------------------------------
CXYMatrix<float> CXYSIS::GetNodeSamples(tree<CXYVector<float> >::iterator  &ritNode, int iSegInd)
{
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = ritNode;
  tree<CXYVector<float> >::iterator root = ptr->begin();
  
  CXYMatrix<float> kMSamplesOrg = GetSamplesOrg();
  kMSamplesOrg = kMSamplesOrg * GetSegLength(iSegInd);
  
  // X, Y, Z, Weight, h1, h2, h3
  float XYZ_n2[3] = {(*pos)[0], (*pos)[1], (*pos)[2]};
  float fWeight = (*pos)[3];
  float fh1 = (*pos)[4];
  float fh2 = (*pos)[5];
  float fh3 = (*pos)[6];
  CXYVector<float> kV_2(3,XYZ_n2);
  
  CXYMatrix<float>& kMSamplesPoints = kMSamplesOrg;
  
  int iNumRow = kMSamplesPoints.GetRows();
  CXYMatrix<float> kMSamplesVector(iNumRow,7);
  for (int iRowInd = 0; iRowInd<iNumRow; iRowInd++) 
  {
    // translation
    CXYVector<float> kVpoint = kMSamplesPoints.GetRow(iRowInd) + kV_2;
    float fVector[7] = {kVpoint[0],kVpoint[1],kVpoint[2],fWeight,fh1,fh2,fh3};
    CXYVector<float> kVector(7,fVector);
    kMSamplesVector.SetRow(iRowInd, kVector);
  }
  
  return kMSamplesVector;
}
//------------------------------------------------------------------------
float CXYSIS::CalBeta_t(CXYMatrix<float> &kM)
{
  int iRow = kM.GetRows();
  CXYVector<float> kV_NextPoint = kM.GetRow(iRow-1);
  
  float h1 = kV_NextPoint[4];
  float h2 = kV_NextPoint[5];
  float fBeta_t =  exp( - (m_fRho_1*h1 + m_fRho_2*h2)/ m_fTau_t );
  
  //  cout << "h1=" << h1 << " ";
  //  cout << "h2=" << h2 << " ";
  //  cout << "Beta_t= " <<  fBeta_t << endl;
  return fBeta_t;
}
//------------------------------------------------------------------------
float CXYSIS::BinSearchConstC(vector<float>& vfBeta_t)
{
  CXYMath<float> kMath;
  sort(vfBeta_t.begin(), vfBeta_t.end());
  float fCons_Min = 0;
  float fCons_Max = 1/vfBeta_t[0];
  
  float fCons_Cur = fCons_Max;
  float fSumofMin = 0;
  
  for (unsigned int i = 0; i<vfBeta_t.size(); i++) 
  {
    //    cout << vfBeta_t[i] << endl;
    fSumofMin = fSumofMin + min(fCons_Cur*vfBeta_t[i],1.0f);
  }
  //  cout << fSumofMin << endl;
  //  cout << "begin" << endl;
  //  cout << fCons_Min << " " << fCons_Cur << " " << fCons_Max << endl;  
  float fCons_Cur_Prev = fCons_Cur;
  
  int iCount = 0;
  while (fabs(fSumofMin-m_iMmax) > kMath.ZERO_TOLERANCE ) 
  {   // for some case it will never stop, so I add condition
    //    cout << fCons_Min << " " << fCons_Cur << " " << fCons_Max << endl;  
    //    if(kMath.AlmostEqual2sComplement(fSumofMin,m_iMmax,10)) 
    //    {
    //        break;
    //    }
    if ((fSumofMin-m_iMmax) > kMath.ZERO_TOLERANCE )
    {
      fCons_Max = fCons_Cur;
    }else{
      fCons_Min = fCons_Cur;
    }
    fCons_Cur = fCons_Min + (fCons_Max - fCons_Min)/2;
    
    fSumofMin = 0;
    for (unsigned int i = 0; i<vfBeta_t.size(); i++) 
    {
      fSumofMin = fSumofMin + min(fCons_Cur*vfBeta_t[i],1.0f);
      //      cout << "ha";
    }
    //    cout << fCons_Cur_Prev << " " << fCons_Cur << endl;
    //    if (fabs(fCons_Cur-fCons_Cur_Prev) < kMath.ZERO_TOLERANCE)
    //    {
    //      break;
    //    }
    if(kMath.AlmostEqual2sComplement(fCons_Cur_Prev,fCons_Cur,1) || iCount > min(max(10,m_iMmax),10))
    {
      break;
    }
    iCount++;
    fCons_Cur_Prev = fCons_Cur;
  }
  cout << "iCount = " << iCount - 1 << endl;
  cout << fSumofMin << endl;
  cout << "vfBeta_t[0]   = " << vfBeta_t[0] << endl;
  cout << "vfBeta_t["<< vfBeta_t.size()<<"]= " << vfBeta_t[vfBeta_t.size()-1] << endl;
  cout << "Cons_Cur=" << fCons_Cur << endl << endl;
  return fCons_Cur;
}
//------------------------------------------------------------------------
float CXYSIS::h1Function(CXYMatrix<float>& kM)
{
  int iNumRow = kM.GetRows();
  int iConInd = iNumRow - 1;

  float fPoint_cur[3] = {kM[iNumRow-1][0],kM[iNumRow-1][1],kM[iNumRow-1][2]};
  CXYVector<float> kVPoint_cur(3,fPoint_cur);
  
  float fCollisionLength = GetCollisionLength();
  
  float fh1 = kM[iNumRow-2][4];
  //  fh1 = 0;
  float fLength = 0;
  float fh1Collision = 0;
  CXYVector<float> vDiff;
  
  vector<CXYTriple<int,int,float> >::iterator tripit;
  CXYVector<float> kVPoint_prv(3);
  for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
    if ((*tripit).second == iConInd && (*tripit).first < (*tripit).second) {

      int ind_prv = (*tripit).first;
      //        float fPoint_prv[3] = {kM[ind_prv][0],kM[ind_prv][1],kM[ind_prv][2]};
      kVPoint_prv[0] = kM[ind_prv][0];
      kVPoint_prv[1] = kM[ind_prv][1];
      kVPoint_prv[2] = kM[ind_prv][2];
      //        CXYVector<float> kVPoint_prv(3,fPoint_prv);
      fLength = (kVPoint_cur - kVPoint_prv).Length();

      fh1Collision = (fLength < fCollisionLength) ? 1 : 0;
      fh1 += fh1Collision;
    }
  }
  
  return fh1;
}

//------------------------------------------------------------------------
float CXYSIS::h2Function(CXYMatrix<float>& kM)
{
  int iNumRow = kM.GetRows();
  int iConInd = iNumRow - 1;
  
  //  cout << iNumRow << endl;
  float fPoint_cur[3] = {kM[iNumRow-1][0],kM[iNumRow-1][1],kM[iNumRow-1][2]};
  CXYVector<float> kVPoint_cur(3,fPoint_cur);
  
  float fh2 = kM[iNumRow-2][5];
  
  //  CXYVector<float> kVContact(iNumRow-1);
  
  float fLength;
  float fh2Prob = 0L;
  
  //int iCount =  GetCountContactNode(iNumRow-1);
  int iCount = GetCountContactNodesFromMap(iNumRow-1);
  
  if (iCount <=3) {
    fh2Prob = 1;
  }else {
    CXYVector<float> kV_Length(iCount);
    CXYVector<float> kV_Connect(iCount);
    
    int iCountInd = 0;
        
    vector<CXYTriple<int,int,float> >::iterator tripit;
    CXYVector<float> kVPoint_prv(3);
    for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
      if ((*tripit).second == iConInd && (*tripit).first < (*tripit).second) {
        int ind_prv = (*tripit).first;
//        float fPoint_prv[3] = {kM[ind_prv][0],kM[ind_prv][1],kM[ind_prv][2]};
        kVPoint_prv[0] = kM[ind_prv][0];
        kVPoint_prv[1] = kM[ind_prv][1];
        kVPoint_prv[2] = kM[ind_prv][2];
//        CXYVector<float> kVPoint_prv(3,fPoint_prv);
        fLength = (kVPoint_cur - kVPoint_prv).Length();
        
        kV_Length[iCountInd] = fLength;
        kV_Connect[iCountInd] = (*tripit).third;
        iCountInd++;
      }
    }
    
    fh2Prob = 1 - (kV_Connect.Dot(kV_Length)/
                   (kV_Connect.Length() * kV_Length.Length()) );
  }
  
  //  cout << kVContact.Length() <<" "<< kVDistInv.Length() << endl;
  //  cout << fh2Prob << endl;
  fh2 = fh2Prob;
  return fh2;
}
//------------------------------------------------------------------------
//int CXYSIS::GetCountContactNode(int iNode)
//{
//  //  CXYMatrix<float> kMObsexp = GetObsexpMatrix();
//  
//  int iRet = 0;
//  for (int i=0; i<iNode; i++) {
//    //    cout << "iNode = " << iNode << " i = " << i << " val = "  << kMObsexp[iNode][i] << endl;
//    if ( (*m_pMDist)[iNode][i] > CXYMath<float>::ZERO_TOLERANCE) {
//      //      cout << kMObsexp[iNode][i] << " ";
//      iRet ++;
//    }
//  }
//  return  (iRet);
//}

void CXYSIS::GetTargetDistribution(void)
{
  tree<CXYVector<float> >* ptr = GetTree();
  int iLevel = ptr->max_depth();
  //  cout << iLevel << endl;
  tree<CXYVector<float> >::fixed_depth_iterator pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  
  pos = ptr->begin_fixed(ptr->begin(),iLevel);
  int iCountSamples = 0;
  while (ptr->is_valid(pos)) 
  {
    iCountSamples ++;
    ++pos;
  }
  //
  
  //new implementation
  pos = ptr->begin_fixed(ptr->begin(),iLevel);
  // enumerate all nodes at same level
  
  int iConformInd = 0;
  CXYVector<float> kV_Weight(iCountSamples);
  CXYVector<float> kV_Collision(iCountSamples);
    
  while (ptr->is_valid(pos)) {
    tree<CXYVector<float> >::iterator prepos = pos;
    CXYMatrix<float> kM(iLevel+1,7);
    int iRowInd = iLevel;
    while (ptr->is_valid(prepos)) {
      kM.SetRow(iRowInd, (*prepos));
      prepos = ptr->parent(prepos);
      iRowInd --;
    }
    
    float fWeight = 0;
//    float fCollision = 0;
    
    //    int iCount = GetCountContactNode(kM.GetRows()-1);
//    int iCount = kM.GetRows();
    
    
    
    
//    cout << iCount << endl;
    
    vector<int>::iterator vit;
    for (vit = m_vConNodeInd.begin(); vit != m_vConNodeInd.end(); vit++) {

      int iNode_1 = (*vit);
      float fPoint_1[3] = {kM[iNode_1][0], kM[iNode_1][1], kM[iNode_1][2]};
      CXYVector<float> kVPoint_1(3,fPoint_1);

      map<int,float> mapConNodeInd = GetContactMap(iNode_1);

      CXYVector<float> kV_Length(mapConNodeInd.size());
      CXYVector<float> kV_Connect(mapConNodeInd.size());

      int iVecInd = 0;

      for (map<int,float>::iterator mapit = mapConNodeInd.begin(); mapit != mapConNodeInd.end(); mapit++) {
        int iNode_2 = (*mapit).first;
        float fPoint_2[3] = {kM[iNode_2][0], kM[iNode_2][1], kM[iNode_2][2]};
        CXYVector<float> kVPoint_2(3,fPoint_2);

        kV_Length[iVecInd] = (kVPoint_2 - kVPoint_1).Length();;
        kV_Connect[iVecInd] = (*mapit).second;
        
//        cout << (*mapit).first <<" " << (*mapit).second << endl;
        iVecInd ++;

      }     
      fWeight = fWeight + 1 - (kV_Connect.Dot(kV_Length)/
                               (kV_Connect.Length() * kV_Length.Length()) );
      
    }
    
    
    
//    for (int iRow=0; iRow<iCount; iRow++) {
//      float fPoint_cur[3] = {kM[iRow][0], kM[iRow][1], kM[iRow][2]};
//      CXYVector<float> kVPoint_cur(3,fPoint_cur);
//      
//      int iCountNodeConnect = GetCountContactNode_All(iRow);
//      
//      CXYVector<float> kV_Length(iCountNodeConnect);
//      CXYVector<float> kV_Connect(iCountNodeConnect);
//      //      cout << iCountNodeConnect << endl;
//      float fLength = 0;
//      int iVecInd = 0;
//      for (int jCol=0; jCol<iCount; jCol++) {
//        if (iRow != jCol && (*m_pMDist)[iRow][jCol] > 0) {
//          float fPoint_prv[3] = {kM[jCol][0], kM[jCol][1], kM[jCol][2]};
//          CXYVector<float> kVPoint_prv(3,fPoint_prv);
//          fLength = (kVPoint_cur - kVPoint_prv).Length();
//          
//          kV_Length[iVecInd] = fLength;
//          kV_Connect[iVecInd] = (*m_pMDist)[iRow][jCol];
////          cout <<iVecInd << "\t" << kV_Length[iVecInd] << "\t" << kV_Connect[iVecInd] << endl;
//          
////          if (fLength < GetCollisionLength() && jCol > iRow) {
////            fCollision ++;
////          }
//          iVecInd ++;
//        }
//      }
//      fWeight = fWeight + 1 - (kV_Connect.Dot(kV_Length)/
//                               (kV_Connect.Length() * kV_Length.Length()) );
////      cout << fWeight << endl;
//    }
////    kV_Weight[iConformInd] = exp(-fWeight/(.1*iCount));
    kV_Weight[iConformInd] = exp(-fWeight);
    
    (*m_pVErr)[iConformInd] = kV_Weight[iConformInd];
    
    iConformInd ++;
    ++pos;
  }
  
//  kV_Weight.Normalize();
  
  pos = ptr->begin_fixed(ptr->begin(),iLevel);
  iConformInd = 0;
  while (ptr->is_valid(pos))
  {
//    (*pos)[3] = (*pos)[3] * kV_Weight[iConformInd];
    (*pos)[3] = kV_Weight[iConformInd];
    ++pos;
    iConformInd++;
  } 
}

//------------------------------------------------------------------------

//int CXYSIS::GetCountContactNode_All(int iNode)
//{
//  //  CXYMatrix<float> kMObsexp = GetObsexpMatrix();
//  int iNumRow = (*m_pMDist).GetRows();
//  int iRet = 0;
//  for (int i=0; i<iNumRow; i++) {
//    if ( (*m_pMDist)[iNode][i] > 0) {
//      iRet ++;
//    }
//  }
//  return  (iRet);
//}

//------------------------------------------------------------------------
void CXYSIS::WritePDB(char* cFN, CXYMatrix<float>& kM, int iInd)
{
  CXYPDB kPDB;
  CXYPDBAtom kAtom;
  
  for (int iRow = 0; iRow < kM.GetRows(); iRow++) 
  {
    kAtom.SetRecName((char *)"ATOM");
    kAtom.SetSerial(iRow);
    kAtom.SetAtomicSym((char *)"CA");
    kAtom.SetResName((char *)"ALA");
    kAtom.SetChainID((char *)"A");
    kAtom.SetResSeq(iRow);
    kAtom.SetX(kM[iRow][0]);
    kAtom.SetY(kM[iRow][1]);
    kAtom.SetZ(kM[iRow][2]);
    kAtom.SetSegID((char *)"C");
    (kPDB.GetPDBAtoms())->push_back(kAtom);
  }
  
  // write comments
  
  // pdb comments 
  string sComment = "";
  char strbuff[256] = "";
  
  sprintf(strbuff,
          "REMARK  100 ErrFunction        = %8.3E\n",
          (*m_pVErr)[iInd]);
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 CollisionFunction  = %8.3E\n",
          (*m_pVEColl)[iInd]);
  sComment = sComment +  strbuff;
  
  sprintf(strbuff, 
          "REMARK  100 PackingDensity     = %8.3f\n",
          GetPackingDensity());
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 CollisionLength    = %8.3f\n", 
          GetCollisionLength());
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 Weight = %8.3E\n",
          kM[kM.GetRows()-1][3]);
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 h1     = %12.3E\n",
          kM[kM.GetRows()-1][4]);
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 h2     = %12.3E\n",
          kM[kM.GetRows()-1][5]);
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 beta_t = %12.3E\n",
          CalBeta_t(kM));
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 Rho1   = %8.3f\n",
          GetRho_1());
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 Rho2   = %8.3f\n",
          GetRho_2());
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 Tau_t  = %8.3f\n",
          GetTau_t());
  sComment = sComment +  strbuff;
  
  sprintf(strbuff,
          "REMARK  100 Adjust = %12.3E\n",
          GetAjust());
  sComment = sComment +  strbuff;
  
  sprintf(strbuff, 
          "REMARK  200 M_max  = %4d\n",
          GetMmax());
  sComment = sComment +  strbuff;
    
  
  
  CXYFile::WriteComment(cFN,"w",sComment.c_str());
  
  // write pdb
  CXYFile::WriteAtomsToPDB(cFN,"a",kPDB);
  
}
//------------------------------------------------------------------------
void CXYSIS::WritePDBArr(void)
{
//  cout << "Writting PDBs..." << endl;
  tree<CXYVector<float> >* ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::fixed_depth_iterator pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  // sort weight
  vector<tree<CXYVector<float> >::iterator> vPosIt;
  vector<pair<float, int> > vWeight;
  int NodeInd = 0;
  float fWeight;
  while (ptr->is_valid(pos)) 
  {
    fWeight = (*pos)[3];
    vWeight.push_back(pair<float,int>(fWeight,NodeInd++));
    vPosIt.push_back(pos);
    ++pos;
  }
  
  //  sort(vWeight.begin(), vWeight.end(), sort_greater_first()); // sort descent
  sort(vWeight.begin(), vWeight.end(), sort_pair_first_greater<float,int>()); // sort descent
  
  char buffer[1024];
  int iScale = 1;
  
  CXYMatrix<float> kM = CXYMatrix<float> (ptr->max_depth() +1,6);
  for (unsigned iit = 0; iit<vWeight.size(); iit++) 
  {
//    cout << vWeight[iit].first << " " << vWeight[iit].second << endl;
    pos = vPosIt[vWeight[iit].second];
    GetOneChain(kM,pos);
    cout << "h1= " << kM[kM.GetRows()-1][4] << endl;
    cout << "h2= " << kM[kM.GetRows()-1][5] << endl;
    cout << "h3= " << kM[kM.GetRows()-1][6] << endl;
    
    for (int iRow=0; iRow < kM.GetRows(); iRow++) 
    {
      for (int iCol=0; iCol < 3; iCol++) 
      {
        kM[iRow][iCol] /= iScale;
      }
    }
    //    cout << "h2= " << kM[kM.GetRows()-1][5] << endl;
    
    //    sprintf(buffer, "%s/%s_%s_chr%s_chr%s_%d_obsexp_%04d.pdb",
    //            GetPDBPath(),"HIC",m_cCell,m_cChrNo,m_cChrNo,m_iResolution,iit+1);
    
   sprintf(buffer, "%s/%04d.pdb",GetOutPath(),iit+1);
   WritePDB(buffer,kM, vWeight[iit].second);
  }
  
}

//------------------------------------------------------------------------
void CXYSIS::GetOneChain(CXYMatrix<float>& rkM,
                                     tree<CXYVector<float > >::iterator itNode)
{
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = itNode;
  
  int iNumRow = ptr->max_depth() + 1;
  
  int iRowInd = iNumRow - 1;
  while (ptr->is_valid(pos)) 
  {
    rkM.SetRow(iRowInd, (*pos));
    --iRowInd;
    pos = ptr->parent(pos);
  }
}
//------------------------------------------------------------------------

void CXYSIS::InitializeChain_SISRoot()
{
  // clear the queue
  m_qPrvNodes.clear();
  
  tree< CXYVector<float> >* ptr;
  tree< CXYVector<float> >::iterator top, root, pos;
  ptr = GetTree();
  
  
  top = ptr->begin();
  float fV0[7] = {0,0,0,1,0,0,0};  // x, y, z, w, h1, h2
  root = ptr->insert(top, CXYVector<float>(7,fV0));
  
  float fLen = GetSegLength(0);
  float fV1[7] = {0,0,fLen,0,0,0,0};
  pos = ptr->append_child(root, CXYVector<float>(7,fV1));
  m_qPrvNodes.push_back(pos);

  //  std::cout << ptr->max_depth() << std::endl;
}

//------------------------------------------------------------------------
void CXYSIS::WritePtsArr(void)
{
//  cout << "Writting Pts..." << endl;
  tree<CXYVector<float> >* ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::fixed_depth_iterator pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  // sort weight
  vector<tree<CXYVector<float> >::iterator> vPosIt;
  vector<pair<float, int> > vWeight;
  int NodeInd = 0;
  float fWeight;
  while (ptr->is_valid(pos)) 
  {
    fWeight = (*pos)[3];
    vWeight.push_back(pair<float,int>(fWeight,NodeInd++));
    vPosIt.push_back(pos);
    ++pos;
  }
  
  //  sort(vWeight.begin(), vWeight.end(), sort_greater_first()); // sort descent
  sort(vWeight.begin(), vWeight.end(), sort_pair_first_greater<float,int>()); // sort descent
  
  char buffer[1024];
  int iScale = 1;
  
  CXYMatrix<float> kM = CXYMatrix<float> (ptr->max_depth() +1,7);
  for (unsigned int iit = 0; iit<vWeight.size(); iit++) 
  {
    cout << vWeight[iit].first << " " << vWeight[iit].second << endl;
    pos = vPosIt[vWeight[iit].second];
    GetOneChain(kM,pos);
    cout << "h1= " << kM[kM.GetRows()-1][4] << endl;
    cout << "h2= " << kM[kM.GetRows()-1][5] << endl;
    cout << "h3= " << kM[kM.GetRows()-1][6] << endl;
    
    for (int iRow=0; iRow < kM.GetRows(); iRow++) 
    {
      for (int iCol=0; iCol < 3; iCol++) 
      {
        kM[iRow][iCol] /= iScale;
      }
    }
    //    cout << "h2= " << kM[kM.GetRows()-1][5] << endl;
    
    //    sprintf(buffer, "%s/%s_%s_chr%s_chr%s_%d_obsexp_%04d.pdb",
    //            GetPDBPath(),"HIC",m_cCell,m_cChrNo,m_cChrNo,m_iResolution,iit+1);
    sprintf(buffer, "%s/%04d.pts",GetOutPath(),iit+1);

    ostream* fp;
    fp = new std::ofstream(buffer);
    for (int iRow = 0; iRow < kM.GetRows(); iRow++) 
    {
      *fp << kM[iRow][0] << "\t" << kM[iRow][1] << "\t" << kM[iRow][2] << "\t"\
      << kM[iRow][3]<< "\t" << kM[iRow][4]<< "\t" << kM[iRow][5] <<"\t" <<kM[iRow][6] << endl;
    }
    delete fp;
    
//    if (iit == 10) {
//      break;
//    }
  }
  
}

//------------------------------------------------------------------------
CXYMatrix<float> CXYSIS::GrowthChain_NoCon(tree< CXYVector<float> >::iterator &ritNode,  CXYVector<float> &rkVPoint)
{
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = ritNode;
  
  int iSize = ptr->depth(pos) + 2;
  CXYMatrix<float> kM(iSize,rkVPoint.GetSize());
  
  // from bottom to top, fill chain 
  int iRowInd = iSize - 1;
  kM.SetRow(iRowInd, rkVPoint);
  
  while (ptr->is_valid(pos)) 
  {
    iRowInd --;
    kM.SetRow(iRowInd, (*pos));
    pos = ptr->parent(pos);
  }
  
  // according previous node h1 and h2, get new h1 and h2
  kM[iSize-1][4] = kM.GetRow(iSize-2)[4];
  kM[iSize-1][5] = kM.GetRow(iSize-2)[5];
  kM[iSize-1][6] = kM.GetRow(iSize-2)[6];
  
  //  char buffer[] = "test.mat";
  //  CXYFile::WriteMatrix("test.mat","w",kM);
  
  return kM;
}

//------------------------------------------------------------------------
// (i, j, value)
// get connection nodes with inode
map<int, float> CXYSIS::GetContactMap(int iNode){
  map<int,float> vMapInd;
  vector<CXYTriple<int,int,float> >::iterator tripit;
  for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
    if ( (*tripit).first == iNode ){
//      cout << " iNode " << iNode << " " << (*tripit).second << endl;
      vMapInd.insert(pair<int,float>(make_pair ((*tripit).second, (*tripit).third)));
    }
    if ( (*tripit).second == iNode ){
//      cout << " iNode " << iNode << " " << (*tripit).first << endl;
      vMapInd.insert(pair<int,float>(make_pair ((*tripit).first, (*tripit).third)));
    }
  }
  return vMapInd;
}

//------------------------------------------------------------------------
int CXYSIS::GetCountContactNodesFromMap(int ind){
  vector<CXYTriple<int,int,float> >::iterator tripit;
  int iNum = 0;
  for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
    if ( (*tripit).second == ind  && (*tripit).first < (*tripit).second ){
      iNum ++;
    }
  }
//  cout << ind << "\t" << iNum << endl;
  return iNum;
}


//------------------------------------------------------------------------
void CXYSIS::RandomPermutation(int ArrN[], int n, int ArrM[], int m){ // m <= n
  MTRand mtrand;
  assert(m<=n);
  for (int i = n-1; i>0; i--) {
    int j = mtrand.randInt(i);
    swap(ArrN[i],ArrN[j]);
  }
  for (int i = 0; i<m; i++) {
    ArrM[i] = ArrN[i];
  }
}

//------------------------------------------------------------------------
CXYVector<float> CXYSIS::GetNodeNextPosition(tree<CXYVector<float> >::iterator  &ritNode, int iSegInd, int ind){
  tree<CXYVector<float> >::iterator pos = ritNode;
  float XYZ_n2[3] = {(*pos)[0], (*pos)[1], (*pos)[2]};
  float fWeight = (*pos)[3];
  float fh1 = (*pos)[4];
  float fh2 = (*pos)[5];
  float fh3 = (*pos)[6];
  CXYVector<float> kV_2(3,XYZ_n2);
  CXYVector<float> kV_NextPoint = (*m_pMSamplesOrg).GetRow(ind) * GetSegLength(iSegInd) + kV_2;
  float fVector[7] = {kV_NextPoint[0],kV_NextPoint[1],kV_NextPoint[2],fWeight,fh1,fh2,fh3};

  CXYVector<float> kV_NextVect(7, fVector );
  return kV_NextVect;
}








//------------------------------------------------------------------------
//CXYMatrix<float> CXYSIS::GrowthChain(tree< CXYVector<float> >::iterator &ritNode, 
//                                     CXYVector<float> &rkVPoint)
//{
//  tree<CXYVector<float> >* ptr = GetTree();
//  tree<CXYVector<float> >::iterator pos = ritNode;
//  
//  int iSize = ptr->depth(pos) + 2;
//  CXYMatrix<float> kM(iSize,6);
//  
//  // from bottom to top, fill chain 
//  int iRowInd = iSize - 1;
//  kM.SetRow(iRowInd, rkVPoint);
//  
//  while (ptr->is_valid(pos)) 
//  {
//    iRowInd --;
//    kM.SetRow(iRowInd, (*pos));
//    pos = ptr->parent(pos);
//  }
//  
//  // according previous node h1 and h2, get new h1 and h2
//  kM[iSize-1][4] = h1Function(kM);
//  kM[iSize-1][5] = h2Function(kM);
//  
//  //  char buffer[] = "test.mat";
//  //  CXYFile::WriteMatrix("test.mat","w",kM);
//  
//  return kM;
//}




//------------------------------------------------------------------------
//------------------------------------------------------------------------

void CXYSIS::GrowOneChain_ByVector(
  tree< CXYVector<float> >::iterator &ritCurNode, 
  tree< CXYVector<float> >::iterator &ritSurfNode,  
  int SegInd, 
  vector<tree<CXYVector<float> >::iterator>* pvChain){

  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = ritCurNode;


  int iSize = SegInd;
//  cout << iSize << endl;
  vector<int> vInd = GetPrvContactsListFromMap(SegInd);
  
//  cout << vInd.size() << endl;

//  copy (vInd.begin(), vInd.end(),
//        ostream_iterator<int>(cout," "));
//  cout << iSize;
//  cout << endl; 
  
  pvChain->push_back(ritSurfNode);

  vector<int>::iterator it;
  for (it = vInd.end()-1; it>=vInd.begin(); it--) {
    while (iSize != *it) {
      iSize --;
      pos = ptr->parent(pos);
    }
    pvChain->push_back(pos);
  }
  
}




//------------------------------------------------------------------------
vector<int> CXYSIS::GetPrvContactsListFromMap(int ind){ 
  vector<CXYTriple<int,int,float> >::iterator tripit;
  vector<int> vInd;
  for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
    if ( (*tripit).second == ind  && (*tripit).first < (*tripit).second ){
      vInd.push_back((*tripit).first);
    }
  }
  // ascending order.
  sort(vInd.begin(), vInd.end());
  return vInd;
}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
CXYVector<float> CXYSIS::GetPrvContactsDistFromMap(int ind){
  vector<CXYTriple<int,int,float> >::iterator tripit;
  vector<pair<int,int> > vIndPair; // (SegInd , index)
  vector<float> vDist;

  int iInd = 0;
  for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
    if ( (*tripit).second == ind  && (*tripit).first < ind ){
      vIndPair.push_back(make_pair((*tripit).first, iInd++));
      vDist.push_back((*tripit).third);
    }
  }
  // ascending order.
  vector<float> vDist_Sorted;
  sort(vIndPair.begin(), vIndPair.end(),RefOrder());
  
  CXYVector<float> kVDist(vDist.size());
  for (unsigned int i=0; i<vIndPair.size(); i++) {
    kVDist[i] = vDist[vIndPair[i].second];
  }
  
  return kVDist;
}

//------------------------------------------------------------------------
void CXYSIS::h1Function_ByVector(vector<tree<CXYVector<float> >::iterator >& rvOneChain){
  
  int iChainConnectedNodeSize = rvOneChain.size();
  CXYVector<float> kV_Vector_cur = *(rvOneChain[0]);
  float fPoint_cur[3] = {kV_Vector_cur[0],kV_Vector_cur[1],kV_Vector_cur[2]};
  CXYVector<float> kV_Point_cur(3,fPoint_cur);
  
  float fCollisionLength_Square = GetCollisionLength() * GetCollisionLength();
  float fLength_Square = 0;
  
  CXYVector<float> kV_Point_prv(3);
  
  float fh1 = (*(rvOneChain[1]))[4];
  for (int i = 1; i < iChainConnectedNodeSize; i ++) {
    
    kV_Point_prv[0] = (*(rvOneChain[i]))[0];
    kV_Point_prv[1] = (*(rvOneChain[i]))[1];
    kV_Point_prv[2] = (*(rvOneChain[i]))[2];
    
    kV_Point_prv -= kV_Point_cur;
    fLength_Square = kV_Point_prv.SquaredLength();
    fh1 += (fLength_Square < fCollisionLength_Square) ? 1 : 0;
  }
  (*(rvOneChain[0]))[4] = fh1;
}

//------------------------------------------------------------------------
void CXYSIS::h2Function_ByVector(vector<tree<CXYVector<float> >::iterator >& rvOneChain, CXYVector<float>& rKV_Dist){
  int iChainConnectedNodeSize = rvOneChain.size();  // n points, 
  int iCount = iChainConnectedNodeSize - 1;        // n-1 edges,
  
  float fh2, fh2Prob ;
  
//  if (iCount <=1) {
//    fh2Prob = 1L;
//  }else {
    CXYVector<float> kV_Vector_cur = *(rvOneChain[0]);
    float fPoint_cur[3] = {kV_Vector_cur[0],kV_Vector_cur[1],kV_Vector_cur[2]};
    CXYVector<float> kV_Point_cur(3,fPoint_cur);

    
    CXYVector<float> kV_Point_prv(3);
    CXYVector<float> kV_Dist_Cal(iCount);
    
    for (int i = iCount,  j=0; i >0; i --, j++) {
      
      kV_Point_prv[0] = (*(rvOneChain[i]))[0];
      kV_Point_prv[1] = (*(rvOneChain[i]))[1];
      kV_Point_prv[2] = (*(rvOneChain[i]))[2];
      
      kV_Point_prv -= kV_Point_cur;
      // ascend sort by connection index
      kV_Dist_Cal[j] = kV_Point_prv.Length();
    }

    fh2Prob = 1 - (kV_Dist_Cal.Dot(rKV_Dist))/(kV_Dist_Cal.Length() * rKV_Dist.Length());
//  }
  fh2 = fh2Prob;
//  cout << "fh2 " << fh2 << endl;
  (*(rvOneChain[0]))[5] = fh2;

}

//------------------------------------------------------------------------
void CXYSIS::h3Function_ByVector(
 vector<tree<CXYVector<float> >::iterator >& rvOneChain, 
 CXYVector<float>& rKV_Dist){
  
  int iChainConnectedNodeSize = rvOneChain.size();  // n points, 
  int iCount = iChainConnectedNodeSize - 1;        // n-1 edges,
  
  float fh3 = 0;
  
  CXYVector<float> kV_Vector_cur = *(rvOneChain[0]);
  float fPoint_cur[3] = {kV_Vector_cur[0],kV_Vector_cur[1],kV_Vector_cur[2]};
  CXYVector<float> kV_Point_cur(3,fPoint_cur);
  
  
  CXYVector<float> kV_Point_prv(3);
  CXYVector<float> kV_Dist_Cal(iCount);
  
  for (int i = iCount,  j=0; i >0; i --, j++) {
    
    kV_Point_prv[0] = (*(rvOneChain[i]))[0];
    kV_Point_prv[1] = (*(rvOneChain[i]))[1];
    kV_Point_prv[2] = (*(rvOneChain[i]))[2];
    
    kV_Point_prv -= kV_Point_cur;
    // ascend sort by connection index
    kV_Dist_Cal[j] = kV_Point_prv.Length();
    if ( abs(kV_Dist_Cal[j] - rKV_Dist[j]) > m_fPersistenceLength ) {
//      cout << iCount << "\t" << kV_Dist_Cal[j] << "\t" << rKV_Dist[j] << endl;
      fh3 ++;
      break;
    }
  }
  
  (*(rvOneChain[0]))[6] = fh3;
  
}
//------------------------------------------------------------------------
void CXYSIS::CalBeta_t_ByVector(vector< vector<tree<CXYVector<float> >::iterator >* >& rvvChains, int iSegInd, vector<float>& rvfBeta_t){
  int iNumChains = rvvChains.size();
  
  // acsend sort by connection index
  CXYVector<float> KV_PrvDist = GetPrvContactsDistFromMap(iSegInd);
  
  float h1, h2, h3;
  float fBeta_t;
  rvfBeta_t.clear();
  for (int iChainInd = 0; iChainInd < iNumChains; iChainInd++) {
      // each chain contain connected nodes position and h1, h2, weight
    vector<tree<CXYVector<float> >::iterator > vOneChain = *(rvvChains[iChainInd]);
    h1Function_ByVector(vOneChain);
    h2Function_ByVector(vOneChain, KV_PrvDist);
    h3Function_ByVector(vOneChain, KV_PrvDist);
    h1 = (*(vOneChain[0]))[4];
    h2 = (*(vOneChain[0]))[5];
    h3 = (*(vOneChain[0]))[6];
    fBeta_t = exp( - (m_fRho_1*h1 + m_fRho_2*h2 + m_fRho_3*h3)/ m_fTau_t );
    rvfBeta_t.push_back(fBeta_t);
  }
}

//------------------------------------------------------------------------
float CXYSIS::GetNoConMaxDist(int iSegInd){
//   vector<CXYTriple<int,int,float> >::iterator tripit;
  // find previous segid and next segid
  int prv_ind = 0;
  int next_ind = 0;
  
  if (iSegInd > (m_vConNodeInd[m_vConNodeInd.size()-1])) {
    prv_ind = m_vConNodeInd[m_vConNodeInd.size()-1];
    next_ind = m_iNumNodes;
  }
  else {
    for (unsigned int i = 0; i < m_vConNodeInd.size(); i++) {
      if (m_vConNodeInd[i] > iSegInd) {
        prv_ind = m_vConNodeInd[i-1];
        next_ind = m_vConNodeInd[i];
        break;
      }
    }
  }
//  cout << prv_ind << " "<< iSegInd << " " <<next_ind << endl;
  // get previous distand and next seg distance
  float prv_dist = 0.0f;
  float next_dist = 0.0f;
  for (int i = prv_ind ; i<iSegInd; i++) {
    prv_dist += GetSegLength(i);
  }
  for (int j = iSegInd; j< next_ind; j++) {
    next_dist += GetSegLength(j);
  }
  
  float lambda = float(2.0f/4.0f);
  prv_dist = prv_dist*(lambda);
  next_dist = next_dist*(lambda);
  
//  cout << min(prv_dist,next_dist) <<endl;

  float nodedist = 0;
  vector<CXYTriple<int,int,float> >::iterator tripit;
  for (tripit = m_vtriple.begin(); tripit != m_vtriple.end(); tripit++) {
    if (((*tripit).first == prv_ind && (*tripit).second == next_ind) ||
        ((*tripit).first == next_ind && (*tripit).second == prv_ind)) {
      nodedist = (*tripit).third;
      break;
    }
  }
  
  float dist_max = 0;
  if (next_ind == m_iNumNodes && m_iNumNodes != m_vConNodeInd[m_vConNodeInd.size()-1]) {
    dist_max = prv_dist+next_dist;
  }else {
    dist_max = min(prv_dist, next_dist)+nodedist;
  }
//  cout << dist_max << " " << nodedist << endl;
  return dist_max;
}

//------------------------------------------------------------------------
CXYVector<float> CXYSIS::GetPrvConPosition(tree<CXYVector<float> >::iterator  &ritNode, int iSegInd){
  // enumerate previous nodes
  tree< CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = ritNode;
  
  int iCount = iSegInd;

  int prv_ind = 0;
  for (unsigned int i = 0; i < m_vConNodeInd.size(); i++) {
    if (m_vConNodeInd[i] > iSegInd) {
      prv_ind = m_vConNodeInd[i-1];
      break;
    }
  }
  if (iSegInd > m_vConNodeInd[m_vConNodeInd.size()-1]) {
    prv_ind = m_vConNodeInd[m_vConNodeInd.size()-1];
  }
  
//  cout << prv_ind << " " << iSegInd << endl;

  while (ptr->is_valid(pos))
  {
    if (iCount == prv_ind) {
//      cout << iCount << " " << (*pos)[0] << " " << (*pos)[1] << " " << (*pos)[2] << endl;
      break;
    }
    pos = ptr->parent(pos);
    iCount --;
  }

  float PrvPoint[] = {(*pos)[0],(*pos)[1],(*pos)[2]};
  CXYVector<float> kvPrvPos = CXYVector<float>(3,PrvPoint);
  return (kvPrvPos);

}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void CXYSIS::SetContIndex(char* cContIndFile){
  CXYMatrix<float> MContInd = CXYFile::ReadMatrix(cContIndFile);
  // Set C++ index
  m_MContInd.SetSize(MContInd.GetRows(),MContInd.GetColumns());
  for (int i=0; i< MContInd.GetRows(); i++) {
    // Segment index
    m_MContInd[i][0] = (int) (MContInd[i][0] - 1);
    // Node start index
    m_MContInd[i][1] = (int) (MContInd[i][1] - 1);
    // Node end index
    m_MContInd[i][2] = (int) (MContInd[i][2] - 1);
//    cout << m_MContInd[i][0] << "\t" 
//          << m_MContInd[i][1] << "\t" 
//          << m_MContInd[i][2] << endl;
  }
}

//------------------------------------------------------------------------
void CXYSIS::SetSegSegPval(char* cPvalFile){
  CXYMatrix<float> MSegSegPval = CXYFile::ReadMatrix(cPvalFile);
  
  // remove i and i+1 connection 
  int iRowCount = 0;
  for (int i=0; i<MSegSegPval.GetRows(); i++) {
    if (MSegSegPval[i][1] - MSegSegPval[i][0] == 1) {
//      cout << "i = " << i <<  ", MSegSegPval[i][0] = " <<  MSegSegPval[i][0] << ", MSegSegPval[i][1] = " <<  MSegSegPval[i][1] << endl  ;
      continue;
    }
    iRowCount ++;
  }
  
  int iColCount = MSegSegPval.GetColumns();
  
  // Set C++ index
  m_MSegSegPval.SetSize(iRowCount,iColCount);
  
  
  for (int i=0; i<iRowCount; i++) {
    if (MSegSegPval[i][1] - MSegSegPval[i][0] == 1) {
      continue;
    }
    // Segment Index
    m_MSegSegPval[i][0] = (int) (MSegSegPval[i][0] - 1);
    // Segment Index
    m_MSegSegPval[i][1] = (int) (MSegSegPval[i][1] - 1);
    // Attraction or Repulsion Index
    // 0 repulsion
    // 2 attraction
    m_MSegSegPval[i][2] = (int) (MSegSegPval[i][2]);
//    cout << m_MSegSegPval[i][0] << "\t" 
//        << m_MSegSegPval[i][1] << "\t" 
//        << m_MSegSegPval[i][2] << endl;
  }
}

//------------------------------------------------------------------------
int CXYSIS::FindSegIndFromNodeInd(int iNodeInd){
  for (int i=0; i<m_MContInd.GetRows(); i++) {
    if (m_MContInd[i][1]<= iNodeInd & iNodeInd <= m_MContInd[i][2]) {
      return m_MContInd[i][0];
    }
  }
  // -1 indicate not found
  return -1;
}

//------------------------------------------------------------------------
int CXYSIS::FindContactionTypeFromSegSegInd(int iSegInd1, int iSegInd2){
  // iSegInd1 less than iSegInd2
  if (iSegInd1 > iSegInd2) {
    int itmp = iSegInd1;
    iSegInd1 = iSegInd2;
    iSegInd2 = itmp;
  }
  
  // 1 4 0 
  // 2 6 2
  // 0 repulsion
  // 2 attraction
  for (int i=0; i<m_MSegSegPval.GetRows(); i++) {
    if (iSegInd1 == m_MSegSegPval[i][0] && iSegInd2 == m_MSegSegPval[i][1]) {
      return m_MSegSegPval[i][2];
    }
  }
  // -1 indicate not found
  return -100;
}

//------------------------------------------------------------------------
vector<int> CXYSIS::GetPrvContactSegIndsFromSegInd(int iSegInd){
  vector<int> vInd;
  for (int i=0; i<m_MSegSegPval.GetRows(); i++) {
    if (m_MSegSegPval[i][1] > iSegInd) {
      break;
    }
    if (m_MSegSegPval[i][1] == iSegInd) {
      vInd.push_back(m_MSegSegPval[i][0]);
    }
  }
  return vInd;
}
//------------------------------------------------------------------------
int CXYSIS::GetPrvContactNumberFromSegInd(int iSegInd){
  vector<int> vInd = GetPrvContactSegIndsFromSegInd(iSegInd);
  return vInd.size();
}


//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_SISAlgorithm(void){

  // Initialization of first segment
  GrowChainSegment_InitRoot();
  
  // Growth for rest segments
//  for (int iSegInd=1; iSegInd<3; iSegInd++) {
  for (int iSegInd=1; iSegInd<m_MContInd.GetRows(); iSegInd++) {
    
//    if (iSegInd == 29) {
//      cout << "" << endl;
//    }
//    if (iSegInd == 10) {
////      cout << "iSegInd = " << iSegInd << endl;
//      break;
//    }
    

    int iCount_SegConnection = GetPrvContactNumberFromSegInd(iSegInd);
    cout << "iSegInd = " << iSegInd <<  " " << iCount_SegConnection << endl;

    int iNodeStartInd, iNodeEndInd;
    FindNodeStartEndIndFromSegInd(iSegInd, iNodeStartInd, iNodeEndInd);
    
    for (int iNodeInd = m_MContInd[iSegInd][1]; iNodeInd<=m_MContInd[iSegInd][2]; iNodeInd++) {

      
      if( iNodeInd % 10 == 0 ) {
        cout << "\tiNodeInd = " << iNodeInd << endl ;  
      }
            
//      if (iCount_SegConnection == 0) {
//        GrowChainSegment_NoRestriction(iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);
//        continue;
//      }
    
      tree<CXYVector<float> > tree_LV;
      tree<CXYVector<float> >::iterator pos_LV;
      // vector to save chains
      vector< vector<tree<CXYVector<float> >::iterator >* > vvChains;
      vector< tree< CXYVector<float> >::iterator > vitPrvNodes;
      GrowChainSegment_GeneratePotentialPoints(
      	tree_LV, pos_LV, vvChains, vitPrvNodes, iNodeInd);

//      cout << "----------------" << iSegInd << " ----------------" << endl;
//      
//      cout << vvChains.size() << endl;
//      for (unsigned int iChainInd =0; iChainInd< vvChains.size(); iChainInd++) {
//        cout << (*( (*(vvChains[iChainInd] )) [0] )) [0] << " " << (*( (*(vvChains[iChainInd] )) [0] )) [1] << " "<<(*( (*(vvChains[iChainInd] )) [0] )) [2] <<  endl ;
//      }
//      cout << "---------------- " << iSegInd << " ----------------" << endl;
//      cout << vitPrvNodes.size() << endl;
//      for (unsigned int iNodeSave =0; iNodeSave < vitPrvNodes.size(); iNodeSave++) {
//        cout <<"Node "<< (*( vitPrvNodes[iNodeSave] ))[0] << " " <<(*( vitPrvNodes[iNodeSave] ))[1]<< " "<<(*( vitPrvNodes[iNodeSave] ))[2] << endl;
//      }
      
    
      // Sequence Importance Sampling
      int iL_t = vvChains.size();
      if (iL_t <= GetMmax()){
        GrowChainSegment_AllSelect (vvChains,vitPrvNodes,iSegInd,iNodeInd,iNodeStartInd,iNodeEndInd);
      }
      else {
        GrowChainSegment_RandSelect(vvChains,vitPrvNodes,iSegInd,iNodeInd,iNodeStartInd,iNodeEndInd);
      }

      // delete the memory we new
      for (unsigned int iNumChainsSize=0; iNumChainsSize <vvChains.size(); iNumChainsSize++) {
        delete vvChains[iNumChainsSize];
      }
      
//      if (iNodeInd == 300) {
//        exit(0);
//      }
      
    } // end for iNodeInd
    

    // cout << ptr->max_depth() << std::endl;
    //    cout << m_MContInd[iSegInd][0] << "\t" << m_MContInd[iSegInd][1] << "\t" << m_MContInd[iSegInd][2] << endl;
  } // end for iSegInd

  //---------------
  // output
  GrowChainSegment_SavePosInd();
  GrowChainSegment_GetErr();
  GrowChainSegment_GetWeight();
  GrowChainSegment_WritePtsArr();

}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_GeneratePotentialPoints(
  tree<CXYVector<float> >& tree_LV,
  tree<CXYVector<float> >::iterator & pos_LV,
  vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
  vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
  int iNodeInd){

  // Nodes to save last points as node
//  tree<CXYVector<float> > tree_LV;
  tree<CXYVector<float> >::iterator root_LV; // , pos_LV;
  CXYVector<float> vRoot(0); // Nothing, just represent root
  tree_LV.set_head(vRoot);
  root_LV = tree_LV.begin();

  tree<CXYVector<float> >::iterator itPrvNode;
  
//  tree< CXYVector<float> >* ptr = GetTree();
  
  
  
  while (!m_qPrvNodes.empty()) {
    itPrvNode = m_qPrvNodes.front();
    m_qPrvNodes.pop_front();
    
    
//    cout << "Prv_x = " << (*itPrvNode)[0] << endl;
//    cout << "Prv_y = " << (*itPrvNode)[1] << endl;
//    cout << "Prv_z = " << (*itPrvNode)[2] << endl;
    
    CXYMatrix<float> kMSamples = GetNodeSamples(itPrvNode,iNodeInd);
    // Conformations for each node
    for (int iSample=0; iSample<kMSamples.GetRows(); iSample++) {
      CXYVector<float> kVPoint = kMSamples.GetRow(iSample);
      // store previous vector, 
      float ftmp[] = {kVPoint[0],kVPoint[1],kVPoint[2], (*itPrvNode)[3] ,(*itPrvNode)[4],(*itPrvNode)[5],(*itPrvNode)[6]};
      CXYVector<float> kV_LV(7,ftmp);
      
//      cout << "PO " << kV_LV[0] <<" " << kV_LV[1] <<" "<<kV_LV[2] <<endl;
      pos_LV = tree_LV.append_child(root_LV, kV_LV);
      
      vector<tree<CXYVector<float> >::iterator > *pvOneChain = new vector<tree<CXYVector<float> >::iterator >;
      //      GrowOneChain_ByVector(itPrvNode, pos_LV, iNodeInd, pvOneChain);
      pvOneChain->push_back(pos_LV);
      vvChains.push_back(pvOneChain);
      // Previous node iterator
      vitPrvNodes.push_back(itPrvNode);
    }
  }
  
//  for (unsigned int iChainInd =0; iChainInd< vvChains.size(); iChainInd++) {
//  	cout << (*( (*(vvChains[iChainInd] )) [0] )) [0] << endl ;
//  }
//  cout << "----------------" << endl;
//  for (unsigned int iNodeSave =0; iNodeSave < vitPrvNodes.size(); iNodeSave++) {
//    cout << (*( vitPrvNodes[iNodeSave] ))[0] << endl;
//  }
  
  
//  cout << vvChains.size() << endl;
}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_RandSelect(
  vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
  vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
  int iSegInd,
  int iNodeInd,
  int iNodeStartInd,
  int iNodeEndInd){
  
  MTRand mtrand;
  tree< CXYVector<float> >* ptr = GetTree();
  tree< CXYVector<float> >::iterator pos;

  int im_t = GetMmax();
  
  // if the node no connection with others, randomly select iMmax node to cont.
  // cout << iNodeInd << " / " << iNumSeg << endl;
  // cout << m_fCollisionLength << endl;
  // calculate beta_t
  
  // vfBeta_t vector include value and index pair
  vector<FloatIntPair> vfBeta_t;

//  CalBeta_t_ByVector(vvChains, iNodeInd, vfBeta_t);
	GrowChainSegment_CalBeta_t(vvChains, vitPrvNodes, iSegInd, iNodeInd, vfBeta_t);

  // binary search const c
  float fConst_C = GrowChainSegment_BinSearchConstC(vfBeta_t);


//  for (vector<FloatIntPair>::iterator it=vfBeta_t.begin(); it!=vfBeta_t.end(); ++it)
//    cout << " " << (*it).first;
//  cout << endl << endl;
//  if (fConst_C == 0) {
//    exit(0);
//  }
//  exit(0);
  
  // get cumsum of min(C*Beta_t,1), order is ascend
  vector<float> vfCumSum, vfMinCBeta1;
  float fCumSum = 0, fMinCBeta1;
  
  for (unsigned int iSample = 0; iSample < vvChains.size(); iSample ++) {
    fMinCBeta1 = min(fConst_C * (vfBeta_t[iSample].first),1.0f);
    fCumSum += fMinCBeta1;
//    cout << "fMinCBeta1 = "<<    fMinCBeta1 << endl;
    vfMinCBeta1.push_back(fMinCBeta1);
    vfCumSum.push_back(fCumSum);
  }
  
  // We can get a real number in the range 0 to 1, excluding 1
  vector<int> vfSampleSelInd;
  float fRand = 0;
  for (int iSameleSel = 1; iSameleSel <= im_t; iSameleSel++) {
    fRand = iSameleSel - mtrand.randExc();
    
    unsigned int iCumSumInd = 0;
    while (fRand >  vfCumSum[iCumSumInd]) {
      if (iCumSumInd == vfCumSum.size()-1) {
        break;
      }
      ++iCumSumInd;
    }
//    cout << "------------------";
//    cout << "fRand = " << fRand << ", vfCumSum[iCumSumInd] = " << vfCumSum[iCumSumInd] << ", iCumSumInd = " << iCumSumInd << ", vfCumSum.size = " << vfCumSum.size();
//    cout << "------------------"<< endl;
    vfSampleSelInd.push_back(vfBeta_t[iCumSumInd].second);
  }
//  cout << "------------------";
//  cout << "vfCumSum[vfCumSum.size()-1] = " << vfCumSum[vfCumSum.size()-1] ;
//  cout << "------------------"<< endl;
  
  // append selected coordinates on certain node
  
  vfSampleSelInd[0] = vfBeta_t[vfBeta_t.size() - 1].second;
//  vfSampleSelInd[1] = vfBeta_t.size() - 2;
//  vfSampleSelInd[2] = vfBeta_t.size() - 3;
//  cout << "vfSampleSelInd[0] = " <<vfSampleSelInd[0] <<  ", vfBeta_t[vfBeta_t.size()-1] = "<< vfBeta_t[vfBeta_t.size()-1] << endl;
  CXYVector<float> kVector(7);
  vector<tree<CXYVector<float> >::iterator > vtmp;
  for (int iSample = 0; iSample < im_t; iSample++) {
    int iChainInd = vfSampleSelInd[iSample];
    kVector = *((*vvChains[iChainInd])[0]);

//    cout << "------------------";
//    cout <<iSegInd << " kV " << iChainInd << " "<< kVector[0] << " " << kVector[1] << " " << kVector[2] <<endl;
//    cout << "------------------";
    
    // recalculate weight
//    cout <<"iChainInd = " << iChainInd<< ", vfMinCBeta1[iChainInd] = " << vfMinCBeta1[iChainInd] << ", kVector[3] = " << kVector[3] <<endl;
    kVector[3] = kVector[3] - log(vfMinCBeta1[iChainInd]);
//    cout <<"kVector[3] = " << kVector[3] <<endl;
    pos = ptr->append_child(vitPrvNodes[iChainInd],kVector);
//cout << "Bt " << kVector[0] << " " << kVector[1] << " " << kVector[2] << endl;
    vtmp.push_back(pos);

    UpdateMap_TreeIt_Vec_Append(vitPrvNodes[iChainInd], pos, iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);

//    UpdateMap_TreeIt_MultiArray_Append(vitPrvNodes[iChainInd], pos, iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);
  }
  UpdateMap_TreeIt_Vec_Clean();
//  UpdateMap_TreeIt_MultiArray_Clean();

  for (int iSample = 0; iSample < im_t; iSample++) {
    m_qPrvNodes.push_back(vtmp[iSample]);
  	m_PrvTreeit_Vec.push_back(vtmp[iSample]);
	}  
}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_InitRoot(void)
{
  // clear the queue
  m_qPrvNodes.clear();
  
  tree< CXYVector<float> >* ptr;
  tree< CXYVector<float> >::iterator top, root, oldpos, newpos;
  ptr = GetTree();
  
  
  top = ptr->begin();
  float fV0[7] = {0,0,0,1,0,0,0};  // x, y, z, w, h1, h2
  root = ptr->insert(top, CXYVector<float>(7,fV0));
  oldpos = root;
  
  float fLen = GetSegLength(0);
  float fV1[7] = {0,0,fLen,1,0,0,0};
  newpos = ptr->append_child(root, CXYVector<float>(7,fV1));

  // Grow first segment
  MTRand mtrand;
  int iNumTrial = max(m_iNumSamplePoints, 100);

  // record previous tree node iterator 
  int iNodeInd;
  for (iNodeInd = 2; iNodeInd<=m_MContInd[0][2]; iNodeInd++) {
  
    oldpos = newpos;

    CXYVector<float> kV_NextVect;
    CXYVector<float> kV_Point(3);
    while (iNumTrial > 0) {
      int iRnd = mtrand.randInt(m_iNumSamplePoints-1);
      kV_NextVect = GetNodeNextPosition(oldpos,iNodeInd,iRnd);
      float Point[] = {kV_NextVect[0],kV_NextVect[1],kV_NextVect[2]};
      kV_Point = CXYVector<float>(3,Point);
      
      if ( ! IsCollision(oldpos, kV_NextVect) ) {
        break;
      }
      
      iNumTrial --;
    }
    
//    if (iNumTrial == 0) {
//      cout << "Give up trial" << endl;
//      exit(0);
//    }

    newpos = ptr->append_child(oldpos,kV_NextVect);
    
  }
  // Only one segment, so need update one time
  // put last point in queue
  m_PrvTreeit_Vec.clear();
  m_PrvTreeit_Vec.push_back(oldpos);
  
//  int iSegInd = FindSegIndFromNodeInd(iNodeInd-1);
  int iSegInd = 0;
  int iNodeStartInd, iNodeEndInd;
  FindNodeStartEndIndFromSegInd(iSegInd, iNodeStartInd, iNodeEndInd);
  

  UpdateMap_TreeIt_Vec_Append(oldpos, newpos, iSegInd, m_MContInd[0][2], iNodeStartInd, iNodeEndInd);
//  UpdateMap_TreeIt_MultiArray_Append(oldpos, newpos, iSegInd, m_MContInd[0][2], iNodeStartInd, iNodeEndInd);

  UpdateMap_TreeIt_Vec_Clean();
//  UpdateMap_TreeIt_MultiArray_Clean();

	m_qPrvNodes.clear();
  m_qPrvNodes.push_back(newpos);
  m_PrvTreeit_Vec.push_back(newpos);
    
 // std::cout << ptr->max_depth() << std::endl;
  
}

//------------------------------------------------------------------------
//void CXYSIS::UpdateOctreeMap_One(
//  tree<CXYVector<float> >::iterator oldpos, 
//  tree<CXYVector<float> >::iterator newpos){
//
//  tree<CXYVector<float> >* ptr = GetTree();
//  tree<CXYVector<float> >::iterator PrvKey = oldpos;
//  
//  Iterator_Int_Map::iterator it = m_Map_PrvKey_Count.find(PrvKey);
//  if (it == m_Map_PrvKey_Count.end()) { // do not find
//    cout << "Error" << endl; 
//    exit(-1);
//  }else { // Count
//    it->second ++;
//    
//    if (m_MapOctree.find(PrvKey) == m_MapOctree.end()) { // if not found, create
//      boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =
//        boost::make_shared< OcTreeNodeAllocator< float , int > >();
//      OcTree<float,int>* pOctree = 
//      new OcTree<float,int>(
//        m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
//      
//      // enumerate the tree to get x, y, z
//      while (ptr->is_valid(oldpos)) 
//      {
//        //cout<< (*oldpos)[0] <<" " <<(*oldpos)[1] <<" " <<(*oldpos)[2] << endl;
//        pOctree->addPoint((*oldpos)[0],(*oldpos)[1],(*oldpos)[2],1,m_maxDepth);
//        oldpos = ptr->parent(oldpos);
//      }
//      m_MapOctree[PrvKey] = pOctree;
//      if (m_MapOctree.find(PrvKey) == m_MapOctree.end()) {
//        cout << "hhhhh" << endl;
//      }
//    }
//    
//    if (it->second == 1) { // First time, append
//      m_MapOctree[newpos] = m_MapOctree[PrvKey];
//      m_MapOctree[newpos]->addPoint(
//        (*newpos)[0],(*newpos)[1],(*newpos)[2],1,m_maxDepth);
//    }else { // Multiple time
//      boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =
//        boost::make_shared< OcTreeNodeAllocator< float , int > >();
//      OcTree<float,int>* pOctree = 
//        new OcTree<float,int>(
//          m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
//      
//      // enumerate the tree to get x, y, z
//      while (ptr->is_valid(oldpos)) 
//      {
//        //cout<< (*oldpos)[0] <<" " <<(*oldpos)[1] <<" " <<(*oldpos)[2] << endl;
//        pOctree->addPoint((*oldpos)[0],(*oldpos)[1],(*oldpos)[2],1,m_maxDepth);
//        oldpos = ptr->parent(oldpos);
//      }
//      pOctree->addPoint((*newpos)[0],(*newpos)[1],(*newpos)[2],1,m_maxDepth);
//      m_MapOctree[newpos] = pOctree;
//    }
//  }
//  
////  cout << m_MapOctree.size() << endl;
//}

//------------------------------------------------------------------------
//void CXYSIS::UpdateOctreeMap_Batch(void){
//  for (Iterator_Int_Map::iterator it=m_Map_PrvKey_Count.begin(); 
//    it!=m_Map_PrvKey_Count.end(); it++) {
//    
//    tree<CXYVector<float> >::iterator PrvKey = it->first;
//    int iCount = it->second;
//
//    if (iCount == 0) {
//      m_MapOctree[PrvKey]->~OcTree();
//      m_MapOctree.erase(PrvKey);
//    } else if (iCount == 1) {
//      m_MapOctree.erase(PrvKey);
//    }
//  }
//  cout << m_MapOctree.size() << endl;
//  m_Map_PrvKey_Count.clear();
//}
//  
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_NoRestriction(int iSegInd, int iNodeInd, int iNodeStartInd, int iNodeEndInd){
  MTRand mtrand;  
  
  vector<tree<CXYVector<float> >::iterator > vtmp;
  tree< CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator oldpos, newpos;
  
  // int iSegInd = FindSegIndFromNodeInd(iNodeInd);
  
  while (!m_qPrvNodes.empty()) {
    oldpos = m_qPrvNodes.front();
    m_qPrvNodes.pop_front(); // pop out
    
    int iNumTrial = max(m_iNumSamplePoints, 100);
    
    CXYVector<float> kV_NextVect;
    CXYVector<float> kV_Point(3);
    while (iNumTrial > 0) {
      int iRnd = mtrand.randInt(m_iNumSamplePoints-1);
      kV_NextVect = GetNodeNextPosition(oldpos,iNodeInd,iRnd);
      float Point[] = {kV_NextVect[0],kV_NextVect[1],kV_NextVect[2]};
      kV_Point = CXYVector<float>(3,Point);
      
      if ( ! IsCollision(oldpos, kV_NextVect) ) {
        break;
      }
      
      iNumTrial --;
    }
    
    newpos = ptr->append_child(oldpos,kV_NextVect);
    vtmp.push_back(newpos);
    // Create new hash key
    UpdateMap_TreeIt_Vec_Append(oldpos, newpos, iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);
//    UpdateMap_TreeIt_MultiArray_Append(oldpos, newpos, iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);

  }
  // Delete old hash key
  UpdateMap_TreeIt_Vec_Clean();
//  UpdateMap_TreeIt_MultiArray_Clean();
  
  // save node to queue
  for (vector<tree<CXYVector<float> >::iterator>::iterator vtmpit = vtmp.begin();
       vtmpit < vtmp.end(); vtmpit++) {
    m_qPrvNodes.push_back(*vtmpit);
    m_PrvTreeit_Vec.push_back(*vtmpit);
  }  
  
}

//------------------------------------------------------------------------
void CXYSIS::UpdateMap_TreeIt_Vec_Append(
                                 tree<CXYVector<float> >::iterator oldpos, 
                                 tree<CXYVector<float> >::iterator newpos,
                                 int iSegInd,
                                 int iNodeInd,
                                 int iNodeStartInd,
                                 int iNodeEndInd){
  
//	if (iSegInd == 0) {return;}
  
  tree<CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator PrvPos = oldpos;

//  cout << endl;
//  while (ptr->is_valid(PrvPos))
//  {
//    cout<< "iSegInd = "<< iSegInd<<", "<< (*PrvPos)[0] <<" " <<(*PrvPos)[1] <<" " <<(*PrvPos)[2] << endl;
//    PrvPos = ptr->parent(PrvPos);
//  }
//  PrvPos = oldpos;
//  cout << endl;

  
  if (m_Treeit_VecpOctree_Map.find(oldpos) == m_Treeit_VecpOctree_Map.end()) {
  
    boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator
    	= boost::make_shared< OcTreeNodeAllocator< float , int > >();
    OcTree<float,int>* pOctree =
	    new OcTree<float,int>(
        m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);

    // add middle points
    CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);
    // Get Middle and End points
    GetMiddlePoints((*PrvPos), (*newpos), kM_MiddlePoints);
    
    for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >=0; iNumPoint --) {
      pOctree->addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
    }
    
    // enumerate the tree to get x, y, z
//    cout << endl;
    while (ptr->is_valid(PrvPos))
    {
//      cout<< (*PrvPos)[0] <<" " <<(*PrvPos)[1] <<" " <<(*PrvPos)[2] << endl;
      pOctree->addPoint((*PrvPos)[0],(*PrvPos)[1],(*PrvPos)[2],1,m_maxDepth);
      PrvPos = ptr->parent(PrvPos);
    }
//    m_MapOctree[PrvKey] = pOctree;
		m_Treeit_VecpOctree_Map[newpos].push_back(pOctree);
  } 
  else {
    

    // first one node of the segment
    m_Treeit_VecpOctree_Map[newpos] = m_Treeit_VecpOctree_Map[oldpos];

//    {
//      if (iSegInd == 3) {
//        float x_ = 0.0f;
//        float y_ = 0.0f;
//        float z_ = 0.0f;
//        float dr = m_fNucleusSphereDiameter;
//        
//        Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
//        Eigen::Matrix<float,4,1> minPosition_ (x_-dr,y_-dr,z_-dr,1);
//        Eigen::Matrix<float,4,1> maxPosition_ (x_+dr,y_+dr,z_+dr,1);
//        vector< OcTreeNode< float, int >* > nodes_;
//        //  prepresentativetarget_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,maxDepth,true);
//        
//        cout << "=================\n";
//        for (int ii = 0; ii< (int) m_Treeit_VecpOctree_Map[newpos].size(); ii++) {
//          m_Treeit_VecpOctree_Map[newpos][ii]->getAllNodesInVolumeOnDepth(
//                                                                                 nodes_,minPosition_,maxPosition_,m_maxDepth,true);
//          
//        }
//        cout << "nodes_.size() = " << nodes_.size() << endl;
//        cout << "m_Treeit_VecpOctree_Map[newpos].size = " << m_Treeit_VecpOctree_Map[newpos].size() << endl;
//        for (unsigned int i=0; i<nodes_.size(); i++) {
//          Eigen::Matrix<float,4,1> neighborpos_ = (nodes_[i])->getPosition();
//          cout << nodes_[i] << "\t";
//          cout << neighborpos_[0] << "," << neighborpos_[1] << "," << neighborpos_[2]  << endl;
//        }
//        
//      }
//    }
                    
    
    OcTree<float,int>* pOctree;
       
    // add middle points including end point
    CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);
    // Get Middle and End points
    GetMiddlePoints((*PrvPos), (*newpos), kM_MiddlePoints);
  
    if (iNodeInd == iNodeStartInd) {
      boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator
    	= boost::make_shared< OcTreeNodeAllocator< float , int > >();

      pOctree = new OcTree<float,int>(
                                     m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
      m_Treeit_VecpOctree_Map[newpos].push_back(pOctree);
      
    } else {

      pOctree = m_Treeit_VecpOctree_Map[newpos][iSegInd];
      
    }

    for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >=0; iNumPoint --) {
      pOctree ->addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
    }
  
  }
  
//  if (iSegInd==10) {
//    cout << "----------------------" << iSegInd << "----------------------\n";
//    for (vector<OcTree<float,int>* >::iterator it= m_Treeit_VecpOctree_Map[newpos].begin(); it != m_Treeit_VecpOctree_Map[newpos].end(); it++) {
//      cout << (*it) << endl;
//    }
//    cout << "----------------------";
//  }
}

void CXYSIS::FindNodeStartEndIndFromSegInd(
                                           int iSegInd, int& iNodeStartInd, int& iNodeEndInd){
  for (int i=0; i<m_MContInd.GetRows(); i++) {
  	if (m_MContInd[i][0] == iSegInd) {
      iNodeStartInd = m_MContInd[i][1];
      iNodeEndInd = m_MContInd[i][2];
      return;
    }
  }
}


void CXYSIS::UpdateMap_TreeIt_Vec_Clean(void){
  for (vector<tree<CXYVector<float> >::iterator>::iterator 
	  it=m_PrvTreeit_Vec.begin(); it != m_PrvTreeit_Vec.end(); it++) {
    TreeIterator_VecpOcTree_Map::iterator tvmit = m_Treeit_VecpOctree_Map.find(*it);
    if ( tvmit != m_Treeit_VecpOctree_Map.end()) {
      m_Treeit_VecpOctree_Map.erase(tvmit);
    }
  }
  m_PrvTreeit_Vec.clear();
}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_AllSelect(
                                        vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                        vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                        int iSegInd,
                                        int iNodeInd,
                                        int iNodeStartInd,
                                        int iNodeEndInd){
  tree< CXYVector<float> >* ptr = GetTree();
  tree< CXYVector<float> >::iterator pos;
                                 
  int im_t = vvChains.size();
  
  vector<tree<CXYVector<float> >::iterator > vtmp;
  for (int iSample=0; iSample<im_t; iSample++) {
    // append child node;
    pos = ptr->append_child(vitPrvNodes[iSample], *( (*vvChains[iSample])[0]));
    vtmp.push_back(pos);
    // update map of treeiterator and vector of poctree
    UpdateMap_TreeIt_Vec_Append(vitPrvNodes[iSample],pos,iSegInd, iNodeInd,iNodeStartInd,iNodeEndInd);
//    UpdateMap_TreeIt_MultiArray_Append(vitPrvNodes[iSample],pos,iSegInd, iNodeInd, iNodeStartInd,iNodeEndInd);
  }
  
  UpdateMap_TreeIt_Vec_Clean();
//  UpdateMap_TreeIt_MultiArray_Clean();
  
  for (int iSample=0; iSample<im_t; iSample++) {
    m_qPrvNodes.push_back(vtmp[iSample]);
    m_PrvTreeit_Vec.push_back(vtmp[iSample]);
  }
  
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_CalBeta_t(
                                        vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                        vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                        int iSegInd,
                                        int iNodeInd,
                                        vector<FloatIntPair>& vBeta_t){
  vBeta_t.clear();

  int iNumChains = (int) vvChains.size();
  int iPrvNodeSegInd = FindSegIndFromNodeInd(iNodeInd-1);
  
  vector<vector<int> > vStartEndInd;
  FindContactStartEndSegIndFromSegInd(iSegInd, vStartEndInd);
  int iAttractionNum = vStartEndInd.size();
  
  float fweight,h1,h2,h3;
  float fBeta_t;
  
  for (int iChainInd =0; iChainInd< iNumChains; iChainInd++) {
    GrowChainSegment_h1Function(
    	*(vvChains[iChainInd]), vitPrvNodes[iChainInd], iSegInd, iNodeInd, iPrvNodeSegInd);
//    cout << "vvv" << iChainInd << " " <<(* (*(vvChains[iChainInd]))[0] ).GetSize() << " " <<
//    	(* (*(vvChains[iChainInd])) [0] )[4]<< endl;
    GrowChainSegment_h2Function(
      *(vvChains[iChainInd]), vitPrvNodes[iChainInd], iSegInd, iNodeInd, iPrvNodeSegInd);
    //    cout << iChainInd << " " <<(* (*(vvChains[iChainInd]))[0] ).GetSize() << " " <<
    //      (* (*(vvChains[iChainInd])) [0] )[5]<< endl;
    if (iAttractionNum != 0) {
//      cout << "iSegInd = " << iSegInd << ", LeftInd = " << vStartEndInd[0][0] << ", RightInd = " << vStartEndInd[0][1] << endl;
      GrowChainSegment_h3Function(
        *(vvChains[iChainInd]), vitPrvNodes[iChainInd], iSegInd, iNodeInd, iPrvNodeSegInd, vStartEndInd);
    }

    fweight = (* (*(vvChains[iChainInd]))[0] )[3];
    h1 = (* (*(vvChains[iChainInd]))[0] )[4] ;
//    h2 = (* (*(vvChains[iChainInd]))[0] )[5] + (*(vitPrvNodes[iChainInd]))[5];
    h2 = (* (*(vvChains[iChainInd]))[0] )[5] ;
    h3 = (* (*(vvChains[iChainInd]))[0] )[6] ;
//    cout << "h1 = " << h1 << ", h2 = " << h2 <<  endl;
    fBeta_t = exp( - (m_fRho_1 * h1 + m_fRho_2 * h2 + m_fRho_3 * h3 )/ m_fTau_t );

//if (iSegInd == 6) {
//  cout << "iSegInd = " << iSegInd << ", iNodeInd = " << iNodeInd << ", fweight = "<< fweight<< ", h1 = " << h1 << ", h2 = " << h2 << ", h3 = " << h3 <<" fBeta_t = " << fBeta_t << endl;
//
//}
    vBeta_t.push_back(FloatIntPair(fBeta_t,iChainInd));
  }

}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_h1Function(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                         tree< CXYVector<float> >::iterator itPrvNode,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iPrvNodeSegInd){
  
	CXYVector<float> kV_Vector_cur = *(vCurChain[0]);
//  cout <<"h1 " << kV_Vector_cur[0] << " " << kV_Vector_cur[1] << " " << kV_Vector_cur[2] <<endl;

//  float fh = (*(vCurChain[0])) [4];
  float fh = 0;
  // 1 is bad, 0 is good
	for (int i=0; i<iSegInd; i++) {
    fh += IsCollision(m_Treeit_VecpOctree_Map[itPrvNode][i],kV_Vector_cur) ? 1 : 0;
  }
//  ( *(vCurChain[0]) )[4] = max(fh,m_Treeit_MultiArray_Map[itPrvNode][iPrvNodeSegInd][1]);
  ( *(vCurChain[0]) )[4] = fh;
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_h2Function(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                         tree< CXYVector<float> >::iterator itPrvNode,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iPrvNodeSegInd){
	CXYVector<float> kV_Vector_cur = *(vCurChain[0]);
//  cout <<iSegInd << " h2 " << kV_Vector_cur[0] << " " << kV_Vector_cur[1] << " " << kV_Vector_cur[2] <<endl;

//  float fh = (*(vCurChain[0])) [5];
  float fh = 0;
  
//  cout << m_Treeit_VecpOctree_Map[itPrvNode].size() << endl;
  vector<int> vSegConInd = GetPrvContactSegIndsFromSegInd(iSegInd);

//  std::ostream_iterator< int > output( cout, " " );
//  std::copy( vSegConInd.begin(), vSegConInd.end(), output );


  // type 0 is repulsion, 2 is attraction, 1 is normal, -1 is not found
  for (vector<int>::iterator it = vSegConInd.begin(); it != vSegConInd.end(); it++ ) {
    int iType = FindContactionTypeFromSegSegInd((*it), iSegInd);
//    cout << iSegInd << " " << *it << " " << iType << endl;
    
    switch (iType) {
      case -100:
        cerr << "not find contact type! " << endl;
        break;
      case -1: // repulsion
        fh += IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][*it], kV_Vector_cur,m_Dr) ? 1 : 0;
//        cout << "Repulsion = " << fh << endl;
        break;
      case 2: // attraction
        fh += IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][*it], kV_Vector_cur,m_Dr) ? 0 : 1;
//          cout << "Attraction = " << fh << endl;
        break;
      default:
        break;
    }
  }
//  cout << "fh = " << fh << endl;

//  if (iSegInd==10) {
//    cout << "----------------------" << iSegInd << "----------------------\n";
//    for (vector<OcTree<float,int>* >::iterator it= m_Treeit_VecpOctree_Map[itPrvNode].begin(); it != m_Treeit_VecpOctree_Map[itPrvNode].end(); it++) {
//      cout << (*it) << endl;
//    }
//    cout << "----------------------";
//  }
  
//  ( *(vCurChain[0]) )[5] = max(fh,m_Treeit_MultiArray_Map[itPrvNode][iPrvNodeSegInd][2]);
  ( *(vCurChain[0]) )[5] = fh;

//  if (fh == 1) {
//    cout << "iSegInd = " << iSegInd<< ", fh = " << fh << endl;
//    // type 0 is repulsion, 2 is attraction, 1 is normal, -1 is not found
//    fh = 0;
//    for (vector<int>::iterator it = vSegConInd.begin(); it != vSegConInd.end(); it++ ) {
//      int iType = FindContactionTypeFromSegSegInd((*it), iSegInd);
//    cout << iSegInd << " " << *it << " " << iType << endl;
//      
//      switch (iType) {
//        case -100:
//          cerr << "not find contact type! " << endl;
//          break;
//        case -1: // repulsion
//                 //        fh += IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][*it], kV_Vector_cur,m_Dr) ? 1 : 0;
//                 //        cout << "Repulsion = " << fh << endl;
//          break;
//        case 2: // attraction
//          fh += IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][*it], kV_Vector_cur,m_Dr) ? 0 : 1;
//          
//          
//          {
////            if (iSegInd == 3) {
//              float x_ = 0.0f;
//              float y_ = 0.0f;
//              float z_ = 0.0f;
//              float dr = m_fNucleusSphereDiameter;
//              
//              Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
//              Eigen::Matrix<float,4,1> minPosition_ (x_-dr,y_-dr,z_-dr,1);
//              Eigen::Matrix<float,4,1> maxPosition_ (x_+dr,y_+dr,z_+dr,1);
//              vector< OcTreeNode< float, int >* > nodes_;
//              //  prepresentativetarget_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,maxDepth,true);
//              
//              cout << "=================\n";
////              for (int ii = 0; ii< (int) m_Treeit_VecpOctree_Map[itPrvNode].size(); ii++) {
//                int ii = *it;
//                m_Treeit_VecpOctree_Map[itPrvNode][ii]->getAllNodesInVolumeOnDepth(
//                         nodes_,minPosition_,maxPosition_,m_maxDepth,true);
//                
////              }
//              cout << "nodes_.size() = " << nodes_.size() << endl;
//              cout << "m_Treeit_VecpOctree_Map[itPrvNode].size = " << m_Treeit_VecpOctree_Map[itPrvNode].size() << endl;
//              for (unsigned int i=0; i<nodes_.size(); i++) {
//                Eigen::Matrix<float,4,1> neighborpos_ = (nodes_[i])->getPosition();
//                cout << nodes_[i] << "\t";
//                cout << neighborpos_[0] << "," << neighborpos_[1] << "," << neighborpos_[2]  << endl;
//              }
//              
////            }
//          }
//          
//          
//          //          cout << "Attraction = " << fh << endl;
//          break;
//        default:
//          break;
//      }
//    }
//  }
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_h3Function(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                         tree< CXYVector<float> >::iterator itPrvNode,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iPrvNodeSegInd,
                                         vector<vector<int > >& vStartEndInd){

	CXYVector<float> kV_Vector_cur = *(vCurChain[0]);
//  cout <<iSegInd << " h3 " << kV_Vector_cur[0] << " " << kV_Vector_cur[1] << " " << kV_Vector_cur[2] <<endl;



//  float fh = ( *(vCurChain[0]) )[6];
  float fh = 0;

  vector<vector<int > >::iterator it, itbegin, itend;
  itbegin = vStartEndInd.begin();
  itend = vStartEndInd.end();
  int iCount, dr;
  int iStartSegInd, iEndSegInd;
  for (it = itbegin; it != itend; it++) {
    iStartSegInd = (*it)[0];
    iEndSegInd = (*it)[1];
    iCount = min(iEndSegInd-iSegInd+1, iSegInd-iStartSegInd);
    dr = iCount * m_Dr;
    fh +=  IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][iStartSegInd], kV_Vector_cur, dr) ? 0 : 1;

  }

  ( *(vCurChain[0]) )[6] = fh;
//  cout << (*(vCurChain[0])).GetSize() << endl;
//  cout << "fh = " << fh << 
//    ", iCount = " << iCount << 
//    ", iSegInd = " << iSegInd <<
//    ", iStartSegInd = " << iStartSegInd <<
//    ", iEndSegInd = " << iEndSegInd << 
//    ", dr = " << dr << endl;
}

//------------------------------------------------------------------------

bool CXYSIS::IsOverlap(OcTree<float,int>* octree_, CXYVector<float>& kV_point, float dr_){
  
  float x_ = kV_point[0];
  float y_ = kV_point[1];
  float z_ = kV_point[2];
  //  cout << "x_ = " << x_ << ", y_ = " << y_ << ", z_ = " << z_ << endl ;
  return IsOverlap(octree_, x_, y_, z_, dr_);
  
}
//------------------------------------------------------------------------
bool CXYSIS::IsOverlap(OcTree<float,int>* pOctree_, float x_, float y_, float z_, float dr_){
  
//  Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
  float min_x_ = max(x_ - dr_, -m_fNucleusSphereDiameter);
  float min_y_ = max(y_ - dr_, -m_fNucleusSphereDiameter);
  float min_z_ = max(z_ - dr_, -m_fNucleusSphereDiameter);
  Eigen::Matrix<float,4,1> minPosition_ (min_x_,min_y_,min_z_,1);
  
  Eigen::Matrix<float,4,1> maxPosition_ (x_+dr_,y_+dr_,z_+dr_,1);
  
  vector< OcTreeNode< float, int >* > nodes_;
  pOctree_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,m_maxDepth,true);

  
  if (nodes_.size() > 0) {
    bool bFlag = false;
    for (int i = 0; i< (int)nodes_.size(); i++) {
      Eigen::Matrix<float,4,1> point_ = (nodes_[i])->getPosition();
      float point_x_ = point_[0];
      float point_y_ = point_[1];
      float point_z_ = point_[2];
      float dist_ = sqrt((point_x_ - x_)*(point_x_ - x_) + (point_y_ - y_ )*(point_y_ - y_) + (point_z_ - z_)*(point_z_ - z_));
      if (dist_ < dr_) {
        bFlag = true;
        break;
        //        cout << x_ << " " << y_ << " " << z_ <<endl;
        //        cout << point_x_ << " " << point_y_ << " " << point_z_ <<endl;
        //        cout << dist_ << endl;
        //        cout << "collision happend";
      }
    }
    return bFlag;  
  }else {
//    cout << "---" << endl;
    return false;
  }
  
}

//------------------------------------------------------------------------
bool CXYSIS::IsCollision(OcTree<float,int>* octree_, CXYVector<float>& kV_point){

  float x_ = kV_point[0];
  float y_ = kV_point[1];
  float z_ = kV_point[2];
//  cout << "x_ = " << x_ << ", y_ = " << y_ << ", z_ = " << z_ << endl ;
  return IsCollision(octree_, x_, y_, z_);
  
}
//------------------------------------------------------------------------
bool CXYSIS::IsCollision(OcTree<float,int>* pOctree_, float x_, float y_, float z_){

//  Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
  float min_x_ = max(x_ - m_dr, -m_fNucleusSphereDiameter);
  float min_y_ = max(y_ - m_dr, -m_fNucleusSphereDiameter);
  float min_z_ = max(z_ - m_dr, -m_fNucleusSphereDiameter);
  
  Eigen::Matrix<float,4,1> minPosition_ (min_x_,min_y_,min_z_,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+m_dr,y_+m_dr,z_+m_dr,1);
  vector< OcTreeNode< float, int >* > nodes_;
  pOctree_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,m_maxDepth,true);
  
  if (nodes_.size() > 0) {
    bool bFlag = false;
    for (int i = 0; i< (int)nodes_.size(); i++) {
      Eigen::Matrix<float,4,1> point_ = (nodes_[i])->getPosition();
      float point_x_ = point_[0];
      float point_y_ = point_[1];
      float point_z_ = point_[2];
      float dist_ = sqrt((point_x_ - x_)*(point_x_ - x_) + (point_y_ - y_ )*(point_y_ - y_) + (point_z_ - z_)*(point_z_ - z_));
      if (dist_ < m_fCollisionLength) {
        bFlag = true;
        break;
        //        cout << x_ << " " << y_ << " " << z_ <<endl;
        //        cout << point_x_ << " " << point_y_ << " " << point_z_ <<endl;
        //        cout << dist_ << endl;
        //        cout << "collision happend";
      }
    }
    return bFlag;  
  }else {
    return false;
  }
  
}
//------------------------------------------------------------------------
float CXYSIS::GrowChainSegment_BinSearchConstC(vector<FloatIntPair>& vfBeta_t)
{
  CXYMath<float> kMath;
  sort(vfBeta_t.begin(), vfBeta_t.end(),FloatIntPairCompare());

  int iSampleSize = (int) vfBeta_t.size();
//  float fInvSampleSize = 1.0f/iSampleSize;

//  if (vfBeta_t[vfBeta_t.size()-1] - vfBeta_t[0] < kMath.ZERO_TOLERANCE ) {
//    cout << "vfBeta_t[0] = " <<vfBeta_t[0] << ", vfBeta_t[vfBeta_t.size()-1]  = "<< vfBeta_t[vfBeta_t.size()-1] << endl;
//    return m_iMmax / (vfBeta_t.size()*vfBeta_t[0]);
//  }

//  cout << "m_iMmax "  << m_iMmax <<endl;
//  cout << "vfBeta_t " << vfBeta_t.size() << endl;
//    for (vector<float>::iterator it=vfBeta_t.begin(); it!=vfBeta_t.end(); ++it)
//      cout << " " << *it;
//  cout << endl;

  float fminnonzero = iSampleSize*m_iMmax;
  int iRest = 1;
  for (int i=0; i<(int)vfBeta_t.size(); i++) {
    if (vfBeta_t[i].first >= kMath.ZERO_TOLERANCE) {
      fminnonzero = vfBeta_t[i].first;
      iRest = (int)vfBeta_t.size()-i;
      break;
    }
  }

  float fCons_Min = 0.0f;
//  float fCons_Max = (float) 1/(vfBeta_t[0]).first;
  float fCons_Max = iRest/fminnonzero;
  
  float ftmp = 0;
  float fSumofMin = (float) (m_iMmax+2);
  float fCons_Cur = fCons_Max;
  float fCons_Prv = fCons_Cur;

  int iMaxCellInd = iSampleSize-1;
  int iPrvCellInd = iMaxCellInd;
  int iCurCellInd = 0;
  
  bool bGreatThanOne = false;
  int iCount = 0;
  vector<float>::iterator it_begin, it_end;

//  while ( fCons_Max - fCons_Cur > fInvSampleSize || 
//         fCons_Cur - fCons_Min > fInvSampleSize) {
//  while (fabs(fSumofMin-m_iMmax)>kMath.ZERO_TOLERANCE && (fCons_Max - fCons_Min)/2 > fInvSampleSize) {
  float iMmax = m_iMmax ;
  while (fabs(fSumofMin-iMmax)>kMath.ZERO_TOLERANCE && (fCons_Max - fCons_Min)/2 > kMath.ZERO_TOLERANCE) {

//  while (iPrvCellInd != iCurCellInd) {
    iPrvCellInd = iCurCellInd;
//    if ( fSumofMin-m_iMmax > CXYMath<float>::ZERO_TOLERANCE )
    if ( fSumofMin-iMmax > 0)
    {
      fCons_Max = fCons_Cur;
    }else{
      fCons_Min = fCons_Cur;
    }
    fCons_Cur = fCons_Min + (fCons_Max - fCons_Min)/2.0f;
    

    fSumofMin = 0;
    bGreatThanOne = false;
    
    for (int i = 0 ; i < iSampleSize ; i++) 
    {
      //    fSumofMin = fSumofMin + min(ftmp,1.0f);
      ftmp = fCons_Cur * vfBeta_t[i].first;
      if (bGreatThanOne) {
        fSumofMin += 1.0f;
      }else {
        if (ftmp < 1.0f) {
          fSumofMin += ftmp;
        }else {
          fSumofMin += 1.0f;
          bGreatThanOne = true;
          iCurCellInd = i;
        }
      }
    }
    
    if (kMath.AlmostEqual2sComplement(fCons_Prv,fCons_Cur,1) || iCount > min(max(100,(int)iMmax), 100) ) {
      break;
    }

    iCount ++;
    fCons_Prv = fCons_Cur;
//    cout << "fCons_Cur = " << fCons_Cur << 
//      ", fSumMin = " << fSumofMin << 
//      ", m_iMmax = " << m_iMmax << 
//      ", MIN = " << fCons_Min << 
//      ", MAX = " << fCons_Max <<endl;
//    cout << "-------------\n";
//    cout << "fSumofMin = " << fSumofMin << ", iMmax =" <<iMmax << ", fCons_Min = " << fCons_Min << ", fCons_Max = " << fCons_Max<< endl;

  }
//  cout << "fSumofMin = " << fSumofMin << ", iMmax =" <<iMmax << ", fCons_Min = " << fCons_Min << ", fCons_Max = " << fCons_Max<< endl;
//  cout << "vfBeta_t[0]   = " << vfBeta_t[0] << endl;
//  cout << "vfBeta_t["<< iSampleSize<<"]= " << vfBeta_t[iSampleSize-1] << endl;
//  cout << "Cons_Cur=" << fCons_Cur << endl << endl;
  return fCons_Cur;
}

//------------------------------------------------------------------------
//void CXYSIS::UpdateMap_TreeIt_MultiArray_Append(tree<CXYVector<float> >::iterator oldpos, 
//                                                tree<CXYVector<float> >::iterator newpos,
//                                                int iSegInd,
//                                                int iNodeInd,
//                                                int iNodeStartInd,
//                                                int iNodeEndInd){
//	if (iNodeInd == 0) {return;}
//    
//  if (m_Treeit_MultiArray_Map.find(oldpos) == m_Treeit_MultiArray_Map.end()) {
//    MultiArray MultiArray_;
////    float fV[] = { (*oldpos)[3], (*oldpos)[4], (*oldpos)[5]};
//    float fV[] = { (*oldpos)[3], (*oldpos)[4], (*oldpos)[5], (*oldpos)[6]};
//    vector<float> v_weight_hvalue = vector<float>(fV,fV+sizeof(fV)/sizeof(float));
//    MultiArray_.push_back(v_weight_hvalue);
//    m_Treeit_MultiArray_Map[newpos] = MultiArray_;
//  } 
//  else {
////    // find previous segment start and end node index
////    int iNodeStartInd, iNodeEndInd;
////    FindNodeStartEndIndFromSegInd(iSegInd, iNodeStartInd, iNodeEndInd);
//    
//    m_Treeit_MultiArray_Map[newpos] = m_Treeit_MultiArray_Map[oldpos];
//    if (iNodeStartInd == iNodeInd ) {
////      float fV[] = { (*oldpos)[3], (*oldpos)[4], (*oldpos)[5]};
////      float fV[] = { (*oldpos)[3], (*oldpos)[4], (*oldpos)[5], (*oldpos)[6]};
//      float fV[] = { (*oldpos)[3], 0.0f, 0.0f, 0.0f};
//      vector<float> v_weight_hvalue = vector<float>(fV,fV+sizeof(fV)/sizeof(float));
//      m_Treeit_MultiArray_Map[newpos].push_back(v_weight_hvalue);
//      // record previous segment last iterator
//      m_PrvSegmentLastTreeit_Vec.push_back(oldpos);
//    }
//    else {
//      // [0] is weight
//      // [1] h1, [2] h2, [3] h3
//      m_Treeit_MultiArray_Map[newpos][iSegInd][0] = m_Treeit_MultiArray_Map[oldpos][iSegInd][0];
//      m_Treeit_MultiArray_Map[newpos][iSegInd][1] = m_Treeit_MultiArray_Map[oldpos][iSegInd][1];
//      m_Treeit_MultiArray_Map[newpos][iSegInd][2] = m_Treeit_MultiArray_Map[oldpos][iSegInd][2];
//      m_Treeit_MultiArray_Map[newpos][iSegInd][3] = m_Treeit_MultiArray_Map[oldpos][iSegInd][3];
//    }
//  }
//}
//------------------------------------------------------------------------
//void CXYSIS::UpdateMap_TreeIt_MultiArray_Clean(){
//
//  vector<tree<CXYVector<float> >::iterator>::iterator it, itbegin, itend;
//  itbegin = m_PrvSegmentLastTreeit_Vec.begin();
//  itend = m_PrvSegmentLastTreeit_Vec.end();
//
//  for ( it = itbegin; it != itend; it++) {
//    TreeIterator_MultiArray_Map::iterator tmmit = m_Treeit_MultiArray_Map.find(*it);
//    if (tmmit != m_Treeit_MultiArray_Map.end()) {
//      m_Treeit_MultiArray_Map.erase(*it);
//    }
//  }
//
//  m_PrvSegmentLastTreeit_Vec.clear();
//
//}
//------------------------------------------------------------------------
void CXYSIS::FindContactStartEndSegIndFromSegInd(int iSegInd, vector<vector<int> > & vStartEndInd){
  int iLeftCloseInd = 0;
  int iRightCloseInd = m_iNumNodes+1 ;
  for (int i=0; i<m_MSegSegPval.GetRows(); i++) {
    if (m_MSegSegPval[i][2] == 2 ) {
      if (m_MSegSegPval[i][0] < iSegInd && m_MSegSegPval[i][1] > iSegInd) {
        iLeftCloseInd = max(iLeftCloseInd, m_MSegSegPval[i][0]);
        iRightCloseInd = min(iRightCloseInd,m_MSegSegPval[i][1]);
      }
    }
  }

  if (iLeftCloseInd !=0 || iRightCloseInd != m_iNumNodes+1) {
    vector<int> vInd;
    vInd.push_back(iLeftCloseInd);
    vInd.push_back(iRightCloseInd);
    vStartEndInd.push_back(vInd);
  }
  
}
//------------------------------------------------------------------------
void CXYSIS::SetContIndex(void){
  m_MContInd.SetSize(m_iNumNodes,3);
  for (int i=0; i<m_iNumNodes; i++) {
    m_MContInd[i][0] = i;
    m_MContInd[i][1] = i;
    m_MContInd[i][2] = i;
  }
}
//------------------------------------------------------------------------
void CXYSIS::GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints){
  
  float dx = (EndPoint[0] - StartPoint[0])/(float)(m_iNumMiddleEndPoints);
  float dy = (EndPoint[1] - StartPoint[1])/(float)(m_iNumMiddleEndPoints);
  float dz = (EndPoint[2] - StartPoint[2])/(float)(m_iNumMiddleEndPoints);
  
  for (int i =0; i< m_iNumMiddleEndPoints; i++) {
    MiddleEndPoints[i][0] = StartPoint[0] + (i+1) * dx;
    MiddleEndPoints[i][1] = StartPoint[1] + (i+1) * dy;
    MiddleEndPoints[i][2] = StartPoint[2] + (i+1) * dz;
  }
  
}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_SavePosInd(void){
  tree<CXYVector<float> > * ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::fixed_depth_iterator fixed_depth_pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  while (ptr->is_valid(fixed_depth_pos)) {
    m_posInd_Vec.push_back(fixed_depth_pos);
    ++fixed_depth_pos;
  }
}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_GetErr(void){
//  cout << "Sorting err value..." << endl;

  tree<CXYVector<float> > * ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::fixed_depth_iterator fixed_depth_pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  tree<CXYVector<float> >::iterator cur_pos, last_pos;
//  int iInd = 0;
  int iChainInd = 0;
  while (ptr->is_valid(fixed_depth_pos)) { // same depth
    last_pos = fixed_depth_pos;
    cur_pos = fixed_depth_pos;
    
    float fh = 0.0f;
    int iSegInd = iLevel -1;
//    cout << "---------------" << endl;
    while (ptr->is_valid(cur_pos) ) { // trace back to root
      vector<int> vSegConInd = GetPrvContactSegIndsFromSegInd(iSegInd);
//      cout << "iSegInd = " << iSegInd << 
//        ", m_Treeit_VecpOctree_Map[last_pos].size() = " << m_Treeit_VecpOctree_Map[last_pos].size() << 
//        ", vSegConInd.size() =" << vSegConInd.size() << endl;

      for (vector<int>::iterator it = vSegConInd.begin(); it != vSegConInd.end(); it++) {
//        cout << "*it = " << *it << " " << m_Treeit_VecpOctree_Map[last_pos][*it]  << endl;
        int iType = FindContactionTypeFromSegSegInd((*it), iSegInd);
        switch (iType) {
          case -1: // repulsion
            fh += IsOverlap(m_Treeit_VecpOctree_Map[last_pos][*it], (*cur_pos),m_Dr) ? 1 : 0;
            //        cout << "Repulsion = " << fh << endl;
            break;
          case 2: // attraction
            fh += IsOverlap(m_Treeit_VecpOctree_Map[last_pos][*it], (*cur_pos),m_Dr) ? 0 : 1;
            //        cout << "Attraction = " << fh << endl;
            break;
          default:
            break;
        }
      }

      iSegInd -- ;
      cur_pos = ptr->parent(cur_pos);
    }
//    cout << fh << " " << iChainInd <<endl;
//    m_Err_Ind.push_back(pair<float,int>(fh,iChainInd));
    m_Err_map[last_pos] = fh;

    iChainInd ++;
    ++fixed_depth_pos;
  }

  // less is good
  sort(m_Err_Ind.begin(), m_Err_Ind.end(),sort_pair_first_less<float,int>()); // asscend
}


void CXYSIS::GrowChainSegment_GetWeight(void){
//  cout << "Sorting weight..." << endl;
  tree<CXYVector<float> >* ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::fixed_depth_iterator pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  int iChainInd = 0;
  while (ptr->is_valid(pos))
  {
    m_Weight_Ind.push_back(pair<float,int>((*pos)[3],iChainInd));
    
    iChainInd++;
    ++pos;
  }
  sort(m_Weight_Ind.begin(), m_Weight_Ind.end(), sort_pair_first_greater<float,int>()); // sort descent
}


void CXYSIS::GrowChainSegment_WritePtsArr(void){
  tree<CXYVector<float> > * ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::iterator pos ;
  
  
  char buffer[1024];
  int iScale = 1;

  
  CXYMatrix<float> kM = CXYMatrix<float> (iLevel +1,7);
//for (int iit=0; iit< (int)m_Err_Ind.size(); iit++) {
//  int iChainInd = m_Err_Ind[iit].second;
  int iNumOutput = min(1, (int)m_Weight_Ind.size());
  iNumOutput = (int)m_Weight_Ind.size();
  for ( int iit = 0; iit<iNumOutput; iit++) 
  {
//    cout << m_Weight_Ind[iit].first << " " << m_Weight_Ind[iit].second << endl;
    int iChainInd = m_Weight_Ind[iit].second;
    pos = m_posInd_Vec[iChainInd];
    GetOneChain(kM,pos);
//    cout << "h1= " << kM[kM.GetRows()-1][4] << " ";
//    cout << "h2= " << kM[kM.GetRows()-1][5] << " ";
//    cout << "h3= " << kM[kM.GetRows()-1][6] << endl;
    
    for (int iRow=0; iRow < kM.GetRows(); iRow++)
    {
      for (int iCol=0; iCol < 3; iCol++) 
      {
        kM[iRow][iCol] /= iScale;
      }
    }
    //    cout << "h2= " << kM[kM.GetRows()-1][5] << endl;
    
    //    sprintf(buffer, "%s/%s_%s_chr%s_chr%s_%d_obsexp_%04d.pdb",
    //            GetPDBPath(),"HIC",m_cCell,m_cChrNo,m_cChrNo,m_iResolution,iit+1);
    sprintf(buffer, "%s/%04d.pts",GetOutPath(),iit+1);
    
//    float fErr = -100.0f;
//    for (int iCount=0; iCount <= m_iMmax; iCount++) {
//      if (m_Err_Ind[iCount].second == iit) {
////        cout << iit << " " << iCount << " " <<m_Err_Ind[iCount].first << endl;
//        fErr = m_Err_Ind[iCount].first;
//        break;
//      }
//    }
    
    ostream* fp;
    fp = new std::ofstream(buffer);
    for (int iRow = 0; iRow < kM.GetRows(); iRow++) 
    {
      *fp << kM[iRow][0] << "\t" << kM[iRow][1] << "\t" << kM[iRow][2] << "\t"\
      << kM[iRow][3]<< "\t" << kM[iRow][4]<< "\t" << kM[iRow][5] <<"\t" <<kM[iRow][6] << endl;
    }
    *fp << "# LogWeight= " << m_Weight_Ind[iit].first << endl;
//    *fp << "# Err= " << m_Err_map[pos] << endl;
    *fp << "# Err= " << m_Err_map[pos] << endl;
    delete fp;
    
  }
  
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_SISAlgorithm_2(void){
  
  // Initialization of first segment
  GrowChainSegment_InitRoot();
  
  // Growth for rest segments
  //  for (int iSegInd=1; iSegInd<3; iSegInd++) {
  for (int iSegInd=1; iSegInd<m_MContInd.GetRows(); iSegInd++) {
    
//    if (iSegInd == 2) {
//      cout << "" << endl;
//    }
    if (iSegInd % 10 == 0) {
      cout << "iSegInd/Total = " << iSegInd << "/" <<m_MContInd.GetRows() <<endl;
    }
//    iSegInd = 6;
    
//    int iCount_SegConnection = GetPrvContactNumberFromSegInd_2(iSegInd);
//    cout << "iSegInd = " << iSegInd <<  " " << iCount_SegConnection << endl;
    
    int iNodeStartInd, iNodeEndInd;
    FindNodeStartEndIndFromSegInd(iSegInd, iNodeStartInd, iNodeEndInd);
//    exit(0);
    
    for (int iNodeInd = m_MContInd[iSegInd][1]; iNodeInd<=m_MContInd[iSegInd][2]; iNodeInd++) {
      
      
//      if( iNodeInd % 10 == 0 ) {
//        cout << "\tiNodeInd = " << iNodeInd << endl ;  
//      }
      
      //      if (iCount_SegConnection == 0) {
      //        GrowChainSegment_NoRestriction(iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);
      //        continue;
      //      }
      
      tree<CXYVector<float> > tree_LV;
      tree<CXYVector<float> >::iterator pos_LV;
      // vector to save chains
      vector< vector<tree<CXYVector<float> >::iterator >* > vvChains;
      vector< tree< CXYVector<float> >::iterator > vitPrvNodes;
      GrowChainSegment_GeneratePotentialPoints(
                                               tree_LV, pos_LV, vvChains, vitPrvNodes, iNodeInd);
      
      //      cout << "----------------" << iSegInd << " ----------------" << endl;
      //      
      //      cout << vvChains.size() << endl;
      //      for (unsigned int iChainInd =0; iChainInd< vvChains.size(); iChainInd++) {
      //        cout << (*( (*(vvChains[iChainInd] )) [0] )) [0] << " " << (*( (*(vvChains[iChainInd] )) [0] )) [1] << " "<<(*( (*(vvChains[iChainInd] )) [0] )) [2] <<  endl ;
      //      }
      //      cout << "---------------- " << iSegInd << " ----------------" << endl;
      //      cout << vitPrvNodes.size() << endl;
      //      for (unsigned int iNodeSave =0; iNodeSave < vitPrvNodes.size(); iNodeSave++) {
      //        cout <<"Node "<< (*( vitPrvNodes[iNodeSave] ))[0] << " " <<(*( vitPrvNodes[iNodeSave] ))[1]<< " "<<(*( vitPrvNodes[iNodeSave] ))[2] << endl;
      //      }
      
      
      // Sequence Importance Sampling
      int iL_t = vvChains.size();
      if (iL_t <= GetMmax()){
        GrowChainSegment_AllSelect (vvChains,vitPrvNodes,iSegInd,iNodeInd,iNodeStartInd,iNodeEndInd);
      }
      else {
        GrowChainSegment_RandSelect_2(vvChains,vitPrvNodes,iSegInd,iNodeInd,iNodeStartInd,iNodeEndInd);
      }
      
      // delete the memory we new
      for (unsigned int iNumChainsSize=0; iNumChainsSize <vvChains.size(); iNumChainsSize++) {
        delete vvChains[iNumChainsSize];
      }
      
      //      if (iNodeInd == 300) {
      //        exit(0);
      //      }
      
    } // end for iNodeInd
    
    
    // cout << ptr->max_depth() << std::endl;
    //    cout << m_MContInd[iSegInd][0] << "\t" << m_MContInd[iSegInd][1] << "\t" << m_MContInd[iSegInd][2] << endl;
  } // end for iSegInd
  
  //---------------
  // output
  GrowChainSegment_SavePosInd();
  GrowChainSegment_GetErr_2();
  GrowChainSegment_GetWeight();
  GrowChainSegment_WritePtsArr();
  
}



//------------------------------------------------------------------------
void CXYSIS::SetSegSegPval_2(char* cPvalFile){
  CXYMatrix<float> MSegSegPval = CXYFile::ReadMatrix(cPvalFile);
  
  // remove i and i+1 connection 
  int iRowCount = 0;
  for (int i=0; i<MSegSegPval.GetRows(); i++) {
    if (MSegSegPval[i][1] - MSegSegPval[i][0] == 1) {
      //      cout << "i = " << i <<  ", MSegSegPval[i][0] = " <<  MSegSegPval[i][0] << ", MSegSegPval[i][1] = " <<  MSegSegPval[i][1] << endl  ;
      continue;
    }
    iRowCount ++;
  }
  
  int iColCount = MSegSegPval.GetColumns();
  
  // Set C++ index
  m_MSegSegPval_2.SetSize(iRowCount,iColCount);
  
  int iCount = 0;
  for (int i=0; i<iRowCount; i++) {
    if ((int)(MSegSegPval[i][1] - MSegSegPval[i][0]) == 1) {
      continue;
    }
    // Segment Index
    m_MSegSegPval_2[iCount][0] = (MSegSegPval[i][0] - 1);
    // Segment Index
    m_MSegSegPval_2[iCount][1] = (MSegSegPval[i][1] - 1);
    // Attraction or Repulsion Index
    // Distance
    m_MSegSegPval_2[iCount][2] = (MSegSegPval[i][2]);
    iCount ++;
//    cout << m_MSegSegPval_2[i][0] << "\t" 
//        << m_MSegSegPval_2[i][1] << "\t" 
//        << m_MSegSegPval_2[i][2] << endl;
  }
}
//------------------------------------------------------------------------
vector<int> CXYSIS::GetPrvContactSegIndsFromSegInd_2(int iSegInd){
  vector<int> vInd;
  for (int i=0; i<m_MSegSegPval_2.GetRows(); i++) {
    if ((int)m_MSegSegPval_2[i][1] > iSegInd) {
      break;
    }
    if ((int)m_MSegSegPval_2[i][1] == iSegInd) {
      vInd.push_back((int)m_MSegSegPval_2[i][0]);
    }
  }
  return vInd;
}

//------------------------------------------------------------------------
int CXYSIS::GetPrvContactNumberFromSegInd_2(int iSegInd){
  vector<int> vInd = GetPrvContactSegIndsFromSegInd_2(iSegInd);
  return vInd.size();
}

//------------------------------------------------------------------------
CXYVector<float> CXYSIS::GetPrvContactsDist_2(int iSegInd){
  vector<float> vInd;
  for (int i=0; i<m_MSegSegPval_2.GetRows(); i++) {
    if ((int)m_MSegSegPval_2[i][1] > iSegInd) {
      break;
    }
    if ((int)m_MSegSegPval_2[i][1] == iSegInd) {
      vInd.push_back(m_MSegSegPval_2[i][2]);
    }
  }
  
  CXYVector<float> kVDist(vInd.size());
  vector<float>::iterator it;
  for (int i=0; i<(int)vInd.size(); i++) {
    kVDist[i] = vInd[i];
  }
  return kVDist;
}
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_RandSelect_2(
                                         vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                         vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iNodeStartInd,
                                         int iNodeEndInd){
  
  MTRand mtrand;
  tree< CXYVector<float> >* ptr = GetTree();
  tree< CXYVector<float> >::iterator pos;
  
  int im_t = GetMmax();
  
  // if the node no connection with others, randomly select iMmax node to cont.
  // cout << iNodeInd << " / " << iNumSeg << endl;
  // cout << m_fCollisionLength << endl;
  // calculate beta_t
  
  // vfBeta_t vector include value and index pair
  vector<FloatIntPair> vfBeta_t;
  
  //  CalBeta_t_ByVector(vvChains, iNodeInd, vfBeta_t);
	GrowChainSegment_CalBeta_t_2(vvChains, vitPrvNodes, iSegInd, iNodeInd, vfBeta_t);


  
  
  
  ////////////////////////////////////////////////////////////
  // binary search const c
  float fConst_C = GrowChainSegment_BinSearchConstC(vfBeta_t);
//  cout << "----------------";
//  cout << fConst_C ;
//  cout << "----------------" << endl;
  
  //  for (vector<FloatIntPair>::iterator it=vfBeta_t.begin(); it!=vfBeta_t.end(); ++it)
  //    cout << " " << (*it).first;
  //  cout << endl << endl;
  //  if (fConst_C == 0) {
  //    exit(0);
  //  }
  //  exit(0);
  
  // get cumsum of min(C*Beta_t,1), order is ascend
  vector<float> vfCumSum, vfMinCBeta1;
  float fCumSum = 0, fMinCBeta1;
  
  for (int iSample = 0; iSample < (int) vvChains.size(); iSample ++) {
    fMinCBeta1 = min(fConst_C * (vfBeta_t[iSample].first),1.0f);
    fCumSum += fMinCBeta1;
    //    cout << "fMinCBeta1 = "<<    fMinCBeta1 << endl;
    vfMinCBeta1.push_back(fMinCBeta1);
    vfCumSum.push_back(fCumSum);
  }
  
  // We can get a real number in the range 0 to 1, excluding 1
  vector<int> vfSampleSelInd;
  vector<float> vfSampleSelMinCBeta1;
  float fRand = 0;
  for (int iSameleSel = 1; iSameleSel <= im_t; iSameleSel++) {
    fRand = iSameleSel - mtrand.randExc();
    
    unsigned int iCumSumInd = 0;
    while (fRand >  vfCumSum[iCumSumInd]) {
      if (iCumSumInd == vfCumSum.size()-1) {
        break;
      }
      ++iCumSumInd;
    }
    //    cout << "------------------";
    //    cout << "fRand = " << fRand << ", vfCumSum[iCumSumInd] = " << vfCumSum[iCumSumInd] << ", iCumSumInd = " << iCumSumInd << ", vfCumSum.size = " << vfCumSum.size();
    //    cout << "------------------"<< endl;
    vfSampleSelInd.push_back(vfBeta_t[iCumSumInd].second);
    vfSampleSelMinCBeta1.push_back(vfMinCBeta1[iCumSumInd]);
  }
  //  cout << "------------------";
  //  cout << "vfCumSum[vfCumSum.size()-1] = " << vfCumSum[vfCumSum.size()-1] ;
  //  cout << "------------------"<< endl;
  ////////////////////////////////////////////////////////////
  
  
  
  // Pick larger beta_t part
//  int iBeta_tInd =vfBeta_t.size() - 1;

  int iBeta_tInd = 0;
  for (int iSampleSel = im_t - 1; 
       iSampleSel >= (int)(im_t*9/10); iSampleSel--, iBeta_tInd++) {
    vfSampleSelMinCBeta1[iSampleSel] = vfMinCBeta1[iBeta_tInd];
    vfSampleSelInd[iSampleSel] = vfBeta_t[iBeta_tInd].second;
  }
  
  // append selected coordinates on certain node
//  vfSampleSelInd[0] = vfBeta_t[vfBeta_t.size() - 1].second;
  //  vfSampleSelInd[1] = vfBeta_t.size() - 2;
  //  vfSampleSelInd[2] = vfBeta_t.size() - 3;
  //  cout << "vfSampleSelInd[0] = " <<vfSampleSelInd[0] <<  ", vfBeta_t[vfBeta_t.size()-1] = "<< vfBeta_t[vfBeta_t.size()-1] << endl;
  CXYVector<float> kVector(7);
  vector<tree<CXYVector<float> >::iterator > vtmp;
  for (int iSample = 0; iSample < im_t; iSample++) {
    int iChainInd = vfSampleSelInd[iSample];
    kVector = *((*vvChains[iChainInd])[0]);
    
    //    cout << "------------------";
    //    cout <<iSegInd << " kV " << iChainInd << " "<< kVector[0] << " " << kVector[1] << " " << kVector[2] <<endl;
    //    cout << "------------------";
    
    // recalculate weight
    //    cout <<"iChainInd = " << iChainInd<< ", vfMinCBeta1[iChainInd] = " << vfMinCBeta1[iChainInd] << ", kVector[3] = " << kVector[3] <<endl;
    kVector[3] = kVector[3] - log(vfSampleSelMinCBeta1[iSample]);
    //    cout <<"kVector[3] = " << kVector[3] <<endl;
    pos = ptr->append_child(vitPrvNodes[iChainInd],kVector);
    //cout << "Bt " << kVector[0] << " " << kVector[1] << " " << kVector[2] << endl;
    vtmp.push_back(pos);
    
    UpdateMap_TreeIt_Vec_Append(vitPrvNodes[iChainInd], pos, iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);
    
    //    UpdateMap_TreeIt_MultiArray_Append(vitPrvNodes[iChainInd], pos, iSegInd, iNodeInd, iNodeStartInd, iNodeEndInd);
  }
  UpdateMap_TreeIt_Vec_Clean();
  //  UpdateMap_TreeIt_MultiArray_Clean();
  
  for (int iSample = 0; iSample < im_t; iSample++) {
    m_qPrvNodes.push_back(vtmp[iSample]);
  	m_PrvTreeit_Vec.push_back(vtmp[iSample]);
	}  
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_CalBeta_t_2(
                                        vector< vector<tree<CXYVector<float> >::iterator >* >& vvChains,
                                        vector< tree< CXYVector<float> >::iterator >& vitPrvNodes,
                                        int iSegInd,
                                        int iNodeInd,
                                        vector<FloatIntPair>& vBeta_t){
  vBeta_t.clear();
  
  int iNumChains = (int) vvChains.size();
  int iPrvNodeSegInd = FindSegIndFromNodeInd(iNodeInd-1);
  
  vector<vector<int> > vStartEndInd;
  FindContactStartEndSegIndFromSegInd_2(iSegInd, vStartEndInd);
  int iAttractionNum = vStartEndInd.size();
  
  vector<int> vConInd = GetPrvContactSegIndsFromSegInd_2(iSegInd);
  CXYVector<float> KV_PrvDist = GetPrvContactsDist_2(iSegInd);

  float fweight,h1,h2,h3;
  float fBeta_t;
  
  for (int iChainInd =0; iChainInd< iNumChains; iChainInd++) {
    GrowChainSegment_h1Function_2(
                                *(vvChains[iChainInd]), vitPrvNodes[iChainInd], iSegInd, iNodeInd, iPrvNodeSegInd);
    //    cout << "vvv" << iChainInd << " " <<(* (*(vvChains[iChainInd]))[0] ).GetSize() << " " <<
    //    	(* (*(vvChains[iChainInd])) [0] )[4]<< endl;
    GrowChainSegment_h2Function_2(
                                *(vvChains[iChainInd]), vitPrvNodes[iChainInd], iSegInd, iNodeInd, iPrvNodeSegInd, vConInd, KV_PrvDist);
    //    cout << iChainInd << " " <<(* (*(vvChains[iChainInd]))[0] ).GetSize() << " " <<
    //      (* (*(vvChains[iChainInd])) [0] )[5]<< endl;
    if (iAttractionNum != 0) {
      //      cout << "iSegInd = " << iSegInd << ", LeftInd = " << vStartEndInd[0][0] << ", RightInd = " << vStartEndInd[0][1] << endl;
      GrowChainSegment_h3Function_2(
                                  *(vvChains[iChainInd]), vitPrvNodes[iChainInd], iSegInd, iNodeInd, iPrvNodeSegInd, vStartEndInd);
    }
    
    fweight = (* (*(vvChains[iChainInd]))[0] )[3];
    h1 = (* (*(vvChains[iChainInd]))[0] )[4] ;
    //    h2 = (* (*(vvChains[iChainInd]))[0] )[5] + (*(vitPrvNodes[iChainInd]))[5];
    h2 = (* (*(vvChains[iChainInd]))[0] )[5] ;
    h3 = (* (*(vvChains[iChainInd]))[0] )[6] ;
    //    cout << "h1 = " << h1 << ", h2 = " << h2 <<  endl;
    fBeta_t = exp( - (m_fRho_1 * h1 + m_fRho_2 * h2 + m_fRho_3 * h3 )/ m_fTau_t );
    
    //if (iSegInd == 6) {
//  cout << "iSegInd= " << iSegInd << " ,iNodeInd= " << iNodeInd << " ,fweight= "<< fweight<< " ,h1= " << h1 << " ,h2= " << h2 << " ,h3 = " << h3 <<" , fBeta_t= " << fBeta_t << endl;

    //}
    vBeta_t.push_back(FloatIntPair(fBeta_t,iChainInd));
  }
  
}

//------------------------------------------------------------------------

void CXYSIS::FindContactStartEndSegIndFromSegInd_2(int iSegInd, vector<vector<int> > & vStartEndInd){
  int iLeftCloseInd = 0;
  int iRightCloseInd = m_iNumNodes+1 ;
  for (int i=0; i<m_MSegSegPval_2.GetRows(); i++) {
    if (m_MSegSegPval_2[i][2] != 0 ) {
      if ((int)m_MSegSegPval_2[i][0] < iSegInd && (int)m_MSegSegPval_2[i][1] > iSegInd) {
        iLeftCloseInd = max(iLeftCloseInd, (int)m_MSegSegPval_2[i][0]);
        iRightCloseInd = min(iRightCloseInd,(int)m_MSegSegPval_2[i][1]);
      }
    }
  }
  
  if (iLeftCloseInd !=0 || iRightCloseInd != m_iNumNodes+1) {
    vector<int> vInd;
    vInd.push_back(iLeftCloseInd);
    vInd.push_back(iRightCloseInd);
    vStartEndInd.push_back(vInd);
  }
  
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_h1Function_2(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                         tree< CXYVector<float> >::iterator itPrvNode,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iPrvNodeSegInd){
  
	CXYVector<float> kV_Vector_cur = *(vCurChain[0]);
  //  cout <<"h1 " << kV_Vector_cur[0] << " " << kV_Vector_cur[1] << " " << kV_Vector_cur[2] <<endl;
  
  //  float fh = (*(vCurChain[0])) [4];
  float fh = 0;
	for (int i=0; i<iSegInd; i++) {
//    fh += IsCollision(m_Treeit_VecpOctree_Map[itPrvNode][i],kV_Vector_cur) ? 1 : 0;
    fh += IsCollision(m_Treeit_VecpOctree_Map[itPrvNode][i],kV_Vector_cur) ? 0 : 1;
  }
  
  if (iSegInd != 0) {
    fh /= iSegInd;
  }
  //  ( *(vCurChain[0]) )[4] = max(fh,m_Treeit_MultiArray_Map[itPrvNode][iPrvNodeSegInd][1]);
  ( *(vCurChain[0]) )[4] = fh;
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_h2Function_2(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                         tree< CXYVector<float> >::iterator itPrvNode,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iPrvNodeSegInd,
                                         vector<int>& rvConInd,
                                         CXYVector<float>& rKV_Dist){
	CXYVector<float> kV_Vector_cur = *(vCurChain[0]);
  tree< CXYVector<float> >* ptr = GetTree();
  tree<CXYVector<float> >::iterator pos = itPrvNode;
  //  cout <<iSegInd << " h2 " << kV_Vector_cur[0] << " " << kV_Vector_cur[1] << " " << kV_Vector_cur[2] <<endl;
  float fPoint_cur[3] = {
    kV_Vector_cur[0], kV_Vector_cur[1], kV_Vector_cur[2]
  };
  CXYVector<float> kV_Point_cur (3,fPoint_cur);
  //  float fh = (*(vCurChain[0])) [5];

  float fh = 0;
  //  cout << m_Treeit_VecpOctree_Map[itPrvNode].size() << endl;
  
  //  std::ostream_iterator< int > output( cout, " " );
  //  std::copy( vSegConInd.begin(), vSegConInd.end(), output );

  int iConNum = rKV_Dist.GetSize();
  if (iConNum != 0) {
    CXYVector<float> kV_Point_prv(3);
    CXYVector<float> kV_Dist_Cal(iConNum);
    
    int iCount = iSegInd -1;
    int iRev = iConNum - 1;
    while (ptr->is_valid(pos)) {
      if (rvConInd[iRev] == iCount) {
//        cout << "iseg = " << iSegInd << ", iCount = " << iCount << endl;
        kV_Point_prv[0] = (*pos)[0];
        kV_Point_prv[1] = (*pos)[1];
        kV_Point_prv[2] = (*pos)[2];
        kV_Point_prv -= kV_Point_cur;
        kV_Dist_Cal[iRev] = kV_Point_prv.Length();  
                      
        iRev --;
        if (iRev == -1) {
          break;
        }
      }
      
      pos = ptr->parent(pos);
      iCount --;
    }

//    float fsigma2 = m_fCollisionLength * m_fCollisionLength;
    //float ffactor = 1/sqrt(2*CXYMath<float>::PI*fsigma2);
    for (int i=0; i<iConNum; i++) {
      float fdiff = abs(kV_Dist_Cal[i] - rKV_Dist[i]);
      //fh += exp( - (fdiff*fdiff)/(2*fsigma2) ) ;
      fh += 1/(1+fdiff);
    }

    fh /= iConNum;
//    if (iConNum==1) {
//      fh = abs(kV_Dist_Cal[0] - rKV_Dist[0]);
//    } else {
//      if (kV_Dist_Cal.Length()==0) { // current point has collision with previous one
//        fh = 1;
//      }else {
//        fh = 1 - (kV_Dist_Cal.Dot(rKV_Dist))/(kV_Dist_Cal.Length()*rKV_Dist.Length());
//      }
//      
//    }
    
    
  }
  
  
//  if (iSegInd==10) {
//    cout << "----------------------" << iSegInd << "----------------------\n";
//    for (vector<OcTree<float,int>* >::iterator it= m_Treeit_VecpOctree_Map[itPrvNode].begin(); it != m_Treeit_VecpOctree_Map[itPrvNode].end(); it++) {
//      cout << (*it) << endl;
//    }
//    cout << "----------------------";
//  }
  
  //  ( *(vCurChain[0]) )[5] = max(fh,m_Treeit_MultiArray_Map[itPrvNode][iPrvNodeSegInd][2]);
  ( *(vCurChain[0]) )[5] = fh;
  
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_h3Function_2(vector<tree<CXYVector<float> >::iterator >& vCurChain,
                                         tree< CXYVector<float> >::iterator itPrvNode,
                                         int iSegInd,
                                         int iNodeInd,
                                         int iPrvNodeSegInd,
                                         vector<vector<int > >& vStartEndInd){
  
	CXYVector<float> kV_Vector_cur = *(vCurChain[0]);
  //  cout <<iSegInd << " h3 " << kV_Vector_cur[0] << " " << kV_Vector_cur[1] << " " << kV_Vector_cur[2] <<endl;
  
  
  
  //  float fh = ( *(vCurChain[0]) )[6];
  float fh = 0;
  
  vector<vector<int > >::iterator it, itbegin, itend;
  itbegin = vStartEndInd.begin();
  itend = vStartEndInd.end();
  int iInsideLoopNum = (int)vStartEndInd.size();
  int iCount, dr;
  int iStartSegInd, iEndSegInd;
  
  
  for (it = itbegin; it != itend; it++) {
    iStartSegInd = (*it)[0];
    iEndSegInd = (*it)[1];
    iCount = min(iEndSegInd-iSegInd+1, iSegInd-iStartSegInd);
    dr = iCount * m_Dr;
//    fh +=  IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][iStartSegInd], kV_Vector_cur, dr) ? 0 : 1;
    fh +=  IsOverlap(m_Treeit_VecpOctree_Map[itPrvNode][iStartSegInd], kV_Vector_cur, dr) ? 1 : 0;
  }
  if (iInsideLoopNum != 0) {
    fh /= iInsideLoopNum;
  }
  
  ( *(vCurChain[0]) )[6] = fh;
  //  cout << (*(vCurChain[0])).GetSize() << endl;
  //  cout << "fh = " << fh << 
  //    ", iCount = " << iCount << 
  //    ", iSegInd = " << iSegInd <<
  //    ", iStartSegInd = " << iStartSegInd <<
  //    ", iEndSegInd = " << iEndSegInd << 
  //    ", dr = " << dr << endl;
}

//------------------------------------------------------------------------
void CXYSIS::GrowChainSegment_GetErr_2(void){
//  cout << "Sorting err value..." << endl;
  
  tree<CXYVector<float> > * ptr = GetTree();
  int iLevel = ptr->max_depth();
  tree<CXYVector<float> >::fixed_depth_iterator fixed_depth_pos = ptr->begin_fixed(ptr->begin(),iLevel);
  
  tree<CXYVector<float> >::iterator cur_pos, last_pos;

  //  int iInd = 0;
  int iChainInd = 0;
  while (ptr->is_valid(fixed_depth_pos)) { // same depth
    last_pos = fixed_depth_pos;
    cur_pos = fixed_depth_pos;
    
    float fh = 0.0f;
    int iSegInd = iLevel -1;
    //    cout << "---------------" << endl;

    while (ptr->is_valid(cur_pos)) {
      float fPoint_cur[3] = {(*cur_pos)[0],(*cur_pos)[1],(*cur_pos)[2]};
      CXYVector<float> kV_Point_cur(3,fPoint_cur);
      
      CXYVector<float> rKV_Dist = GetPrvContactsDist_2(iSegInd);
      vector<int> rvConInd = GetPrvContactSegIndsFromSegInd_2(iSegInd);
      int iConNum = rKV_Dist.GetSize();


      if (iConNum != 0) {
        CXYVector<float> kV_Point_prv(3);
        CXYVector<float> kV_Dist_Cal(iConNum);
        
        int iCount = iSegInd -1;
        int iRev = iConNum - 1;
        tree<CXYVector<float> >::iterator prv_pos = ptr->parent(cur_pos);
        
        while (ptr->is_valid(prv_pos)) {
          if (rvConInd[iRev] == iCount) {
            //        cout << "iseg = " << iSegInd << ", iCount = " << iCount << endl;
            kV_Point_prv[0] = (*prv_pos)[0];
            kV_Point_prv[1] = (*prv_pos)[1];
            kV_Point_prv[2] = (*prv_pos)[2];
            kV_Point_prv -= kV_Point_cur;
            kV_Dist_Cal[iRev] = kV_Point_prv.Length();  
            
            
            iRev --;
            if (iRev == -1) {
              break;
            }
          }
          
          prv_pos = ptr->parent(prv_pos);
          iCount --;
        }

        //float ffactor = 1/sqrt(2*CXYMath<float>::PI*fsigma2);
        for (int i=0; i<iConNum; i++) {
          float fdiff = abs(kV_Dist_Cal[i] - rKV_Dist[i]);
          fh += fdiff;
        }
        
//        fh /= iConNum;
        
//        if (iConNum==1) {
//          fh = abs(kV_Dist_Cal[0] - rKV_Dist[0]);
//        } else {
//          
//          if (kV_Dist_Cal.Length() == 0) { // current point has collision with previous one
//            fh += 1;
//          }else {
//            fh += 1 - (kV_Dist_Cal.Dot(rKV_Dist))/(kV_Dist_Cal.Length()*rKV_Dist.Length());
//          }
//        }

//        cout << "iSegInd = "<<iSegInd << ", fh =" << fh<<endl;
      }


      iSegInd --;
      if (iSegInd == 0) {
        break;
      }
      cur_pos = ptr->parent(cur_pos);
    }
    
    //    cout << fh << " " << iChainInd <<endl;
    //    m_Err_Ind.push_back(pair<float,int>(fh,iChainInd));
    m_Err_map[last_pos] = fh;
    
    iChainInd ++;
    ++fixed_depth_pos;
  }
  
  // less is good
  sort(m_Err_Ind.begin(), m_Err_Ind.end(),sort_pair_first_less<float,int>()); // asscend
}
