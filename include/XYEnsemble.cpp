#include "XYEnsemble.h"
#include "XYSIS.h"
#include <vector>
#include <time.h>       /* time */
//#include "debug_new.h"

using namespace spatialaggregate;
using std::vector;

//----------------------------------------------------------------------------
CXYEnsemble::CXYEnsemble()
{
}

CXYEnsemble::CXYEnsemble(
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
    )
{
	m_cOutPath = new char[1024];
	SetOutPath(cOutPath);
	SetPersistenceLength( fPersistenLength );
	SetPackingDensity(fPackingDensity);
	SetNucleusSphereDiameter(fNucleusSphereDiameter);
	SetNumChains(NumberofChains);
	SetChainLengths(ChainLengths);

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
		SetSegLengths(cStartEndFile, "AVG");
		SetNumNodes(m_vfSegLength.size()); // set node numbers
    SetContIndex(cContIndFile);
	}
	SetCollisionLength( fCollisionLength );

    mNum =1;   //number of the chains

	//---------------------------------------- GAMZE 10/2 -------------------------------------




        for (int i = 0; i < mNum; i++)
        {

            vector <tree <CXYVector <float> >* >  v;

            for (int j=0; j<(2*m_fNumofChains); j++)
            {

                tree <CXYVector<float> >* TreeTemp= new tree< CXYVector <float> >  () ;
                v.push_back(TreeTemp);
                alpha[j] = 0;
                for (int k=0; k<3; k++)
                {
                    endpoints[i][j][k]=0;
                    centromeres[i][j][k]=0;
                }
            // delete TreeTemp;
            }

         m_pTChain.push_back(v);
    }







	//--------------------------------------------------------------------------------------------------

	m_pMSamplesOrg = new CXYMatrix<float>;
	SetSamplesOrg();

  m_iNumMiddleEndPoints = max(int(m_fPersistenceLength / m_fCollisionLength)-1 ,1);

  cout << "persistence length = " << m_fPersistenceLength << endl;
  cout << "collision length = " << m_fCollisionLength << endl;
  cout << "number nodes = " << m_iNumNodes << endl;
  cout << "m_iNumMiddleEndPoints = " << m_iNumMiddleEndPoints <<endl;


  // Construct octree
  m_center = Eigen::Matrix< float, 4, 1 > (0.0f,0.0f,0.0f,1);
  m_minimumVolumeSize = m_fCollisionLength/2;
  m_dr = m_fCollisionLength*2;
//  m_dr = m_fNucleusSphereDiameter/2;

//----------------------------------------------GAMZE 10/2--------------------------------------------------------------------------------------------------------------

    for (int i = 0; i < mNum; i++)
	{
        boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
        m_pOctree.push_back(new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator) );

	}

	//logweight per sample

	  for (int i=0 ; i<mNum; i++)
        {
            //cout<<"index is "<<i<<endl;
            ms_logweight.push_back(FloatIntPair(0,i));

            //cout << "logweight is"<<m_LogWeight[i].first <<endl;
        }

        //logweight per chain

//  for (int i=0 ; i<m_fNumofChains; i++)
//        {
//            //cout<<"index is "<<i<<endl;
//            m_LogWeight.push_back(FloatIntPair(0,i));
//
//            //cout << "logweight is"<<m_LogWeight[i].first <<endl;
//        }

m_maxDepth = ceil(m_pOctree[0]->depthForVolumeSize(m_minimumVolumeSize));
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cout << "minimumVolumeSize =" << m_minimumVolumeSize << endl;
  cout << "maxDistance =" << m_fNucleusSphereDiameter << endl;
  cout << "maxDepth = " <<  m_maxDepth << endl;


}

//----------------------------------------------------------------------------
CXYEnsemble::~CXYEnsemble(void)
{


    for (int i=0 ; i<mNum; i++)
    {
        delete  m_pOctree[i];
   //     delete  m_pTChain1[i];

    }



 for (int i=0 ; i<mNum; i++)
    {
        for (int j=0; j<(2*m_fNumofChains); j++)
        {
            delete m_pTChain[i][j];

        }
            //delete  v[i];

    }


	delete m_cOutPath;
	delete m_pMSamplesOrg;
	//m_pOctree.empty();
    //delete  m_fChainLengths;//.empty();
//    m_pTChain.empty();
  //  m_LogWeight.empty();
  //  ms_logweight.empty();
}

//----------------------------------------------------------------------------
void CXYEnsemble::SetPersistenceLength(float fPL)
{
	m_fPersistenceLength = fPL;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetPersistenceLength(void)
{
	return m_fPersistenceLength;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNumChains(int NofC)
{
	m_fNumofChains = NofC;
}
//----------------------------------------------------------------------------
int CXYEnsemble::GetNumChains(void)
{
	return m_fNumofChains;
}
//------------------------------------------------------------------------------

void CXYEnsemble::SetChainLengths(vector <int> fCLs)
{


    int NofC = GetNumChains();

    CXYVector <int> cent = getCentromeres(NofC);

    int ind = 0;

    for (int i=0; i<2*NofC; i++)
    {
        m_fChainLengths.push_back(0);
    }

    for (int i=0; i<NofC; i++)
    {
        m_fChainLengths[ind]=cent[i];
        m_fChainLengths[ind+1]=(fCLs[i]-cent[i]+1);

        if (ind == 22)
        {
            m_fChainLengths[ind+1]=fCLs[i]-28;
            m_fChainLengths[ind]=29;
        }

        ind = ind +2;
    }

}
//----------------------------------------------------------------------------
vector <int> CXYEnsemble::GetChainLengths(void)
{

     int NofC = GetNumChains();

	return m_fChainLengths;
}
//------------------------------------------------------------------------------
void CXYEnsemble::SetCollisionLength(float fCollision)
{
	vector<float>& vfSegment = GetSegLengths();
	float min_diameter = *(min_element(vfSegment.begin(), vfSegment.end()));
	m_fCollisionLength = min(fCollision,min_diameter);
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetCollisionLength(void)
{
	return m_fCollisionLength;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetPackingDensity(float fPackingDensity)
{
	m_fPackingDensity = fPackingDensity;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetPackingDensity(void)
{
	return m_fPackingDensity;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNucleusSphereDiameter(float fNucleusSphereDiameter)
{
	m_fNucleusSphereDiameter = fNucleusSphereDiameter;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetNucleusSphereDiameter(void)
{
	return m_fNucleusSphereDiameter;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetSBPSphereDiameter(void)
{
    float SBPdiameter = 6000;
	return SBPdiameter;
}
//------------------------------------------------------------------------------
void CXYEnsemble::SetNumNodes(int iNumNodes){
	m_iNumNodes = iNumNodes;
}
//----------------------------------------------------------------------------
int CXYEnsemble::GetNumNodes(void){
	return m_iNumNodes;
}

//----------------------------------------------------------------------------
// Generate sphere sample points with radius persistence length.
void CXYEnsemble::SetSamplesOrg(void)
{
	CXYSO3Sequence sO3sequence(m_iNumSamplePoints);
	sO3sequence.SetSO3Sequence();
	(*m_pMSamplesOrg) = sO3sequence.GetSO3Sequence();
}
//----------------------------------------------------------------------------
CXYMatrix<float> CXYEnsemble::GetSamplesOrg()
{
	return (*m_pMSamplesOrg);
}

//----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::NewXYZ(CXYVector<float> kV_0, CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//	CXYMatrix<float> kMXYZ(3,3);
//	CXYVector<float> kV_Z = kV_2 - kV_1;
//	CXYVector<float> kV_01 = kV_1 - kV_0;
//	kV_01.Normalize();
//	if (fabs(kV_Z.Dot(kV_01)) < CXYMath<float>::ZERO_TOLERANCE)
//	{ // two line parallel
//		kMXYZ = NewXYZ(kV_1,kV_2);
//	} else {
//		float a = kV_Z[0],  b = kV_Z[1],  c = kV_Z[2];
//		float x = kV_01[0], y = kV_01[1], z = kV_01[2];
//
//		float fX[3] = {	y*c - z*b, z*a - x*c, x*b - y*a };
//		CXYVector<float> kV_X(3,fX);
//		kV_X.Normalize();
//
//		x = kV_Z[0], y = kV_Z[1], z = kV_Z[2];
//		a = kV_X[0], b = kV_X[1], c = kV_X[2];
//		float fY[3] = {	y*c - z*b, z*a - x*c, x*b - y*a };
//		CXYVector<float> kV_Y(3,fY);
//		kV_Y.Normalize();
//		kMXYZ.SetRow(0, kV_X);
//		kMXYZ.SetRow(1, kV_Y);
//		kMXYZ.SetRow(2, kV_Z);
//	}
//
//	return kMXYZ;
//}
////----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::NewXYZ(CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//	CXYMatrix<float> kMXYZ(3,3);
//	CXYVector<float> kV_Z = kV_2 - kV_1;
//	kV_Z.Normalize();
//	float a = kV_Z[0], b= kV_Z[1], c = kV_Z[2];
//
//	float x, y, z;
//	if (c > CXYMath<float>::ZERO_TOLERANCE )
//	{
//		x = 1; y = 0; z = -(a*x+b*y)/c;
//	} else if (b > CXYMath<float>::ZERO_TOLERANCE )
//	{
//		x = 0; z = 1; y = -(a*x+c*z)/b;
//	} else {
//		y = 1; z = 0; x = -(b*y+c*z)/a;
//	}
//	float fX[3] = {x,y,z};
//	CXYVector<float> kV_X(3,fX);
//	float fY[3] = {b*z-c*y, c*x-a*z, a*y-b*x};
//	CXYVector<float> kV_Y(3,fY);
//
//	kMXYZ.SetRow(0, kV_X);
//	kMXYZ.SetRow(1, kV_Y);
//	kMXYZ.SetRow(2, kV_Z);
//
//
//	return kMXYZ;
//}

//----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::GetRotMatrix(CXYMatrix<float> &rkMXYZ)
//{
//	CXYVector<float> kV_X = rkMXYZ.GetRow(0);
//	CXYVector<float> kV_Y = rkMXYZ.GetRow(1);
//	CXYVector<float> kV_Z = rkMXYZ.GetRow(2);
//
//	CXYMatrix<float> kM_A(3,3);
//	float alpha, beta, gamma;
//	float Z3diff = fabs(1-kV_Z[2]*kV_Z[2]);
//	if ( Z3diff < CXYMath<float>::ZERO_TOLERANCE)
//	{
//		alpha = acos( min(float(1.0), max(float(-1.0), kV_X[0])));
//		beta  = 0;
//		gamma = 0;
//	} else {		// http://www.macosxguru.net/article.php?story=20040210124637626
//		float Denorm = sqrt(Z3diff);
//		alpha = acos(min(float(1.0), max(float(-1.0), -kV_Z[1]/Denorm)));
//		beta  = acos(min(float(1.0), max(float(-1.0), -kV_Z[2])));
//		gamma = acos(min(float(1.0), max(float(-1.0), -kV_Y[2]/Denorm)));
//	}
//
//	kM_A[0][0] = cos(gamma) * cos(alpha) - cos(beta) * sin(alpha) * sin(gamma);
//	kM_A[0][1] = cos(gamma) * sin(alpha) + cos(beta) * cos(alpha) * sin(gamma);
//	kM_A[0][2] = sin(gamma) * sin(beta);
//	kM_A[1][0] = -sin(gamma) * cos(alpha) - cos(beta) * sin(alpha) * cos(gamma);
//	kM_A[1][1] = -sin(gamma) * sin(alpha) + cos(beta) * cos(alpha) * cos(gamma);
//	kM_A[1][2] = cos(gamma) * sin(beta);
//	kM_A[2][0] = sin(beta) * sin(alpha);
//	kM_A[2][1] = -sin(beta) * cos(alpha);
//	kM_A[2][2] = cos(beta);
//
//	kM_A.GetInverse(kM_A);
//
//
//	return kM_A;
//}

CXYVector <int> CXYEnsemble::getCentromeres(int numofchr)
{
    CXYVector <int> centr(numofchr);

    //enter the mononers that have centromeres

    centr[0]=10;
    centr[1]=15;
    centr[2]=7;
    centr[3]=29;
    centr[4]=10;
    centr[5]=9;
    centr[6]=33;
    centr[7]=7;
    centr[8]=23;
    centr[9]=29;
    centr[10]=29;
    centr[11]=10;
    centr[12]=17;
    centr[13]=41;
    centr[14]=21;
    centr[15]=21;

    return centr;
}

//----------------------------------------------------------------------------
void CXYEnsemble::InitializeChain(int m_samples)
{




  //---------------------------------------------------GAMZE 10/2-----------------------------------------------

//        for ( int i = 0; i < m; i++)
//        {
//           // cout<<"index is "<<i<<endl;
//            m_LogWeight[i].first=0;
//
//           // cout << "logweight is"<<m_LogWeight[i].first <<endl;
//        }

        for ( int i = 0; i < m_samples; i++)
        {
          //  cout<<"index is "<<i<<endl;
            ms_logweight[i].first=0;
            ms_logweight[i].second = i;

          //  cout << "logweight is"<<m_LogWeight[i].first <<endl;
        }


//if octree is not null, then delete it ------------------
  for ( int i = 0; i < m_pOctree.size(); i++ )
  {

       if (m_pOctree[i] -> root_) {
            delete m_pOctree[i];

       }
  }
m_pOctree.clear();
  //---------------------------------------------------------------------------------------------------------------------


  OcTree<float,int> * p_octree;



//----------------------here I am-------------------------


 // vector < vector < tree < CXYVector <float> >*> > ptr_right = GetTree_right();
  vector < vector < tree < CXYVector <float> >*> > ptr = GetTree();



	// erase all nodes if not empty


//
//for (int k=0; k<2; k++)
//{
for (int j=0; j < m_samples; j++)
{
  for (int i = 0; i < (2*m_fNumofChains); i++)
  {
      if ( !ptr[j][i]-> empty() )
      {
              //  cout<<"initialize chain"<<endl;

          ptr[j][i]-> clear();
         // cout<<"burda "<<endl;

      }
  }
}
//}


//vector < vector  <  tree < CXYVector <float> >::iterator> >  top_r, root_r;
vector < vector  <  tree < CXYVector <float> >::iterator> >  top, root;



for (int j=0; j<m_samples; j++)
{
    vector  <  tree < CXYVector <float> >::iterator> t1,r1;

    for (int i = 0; i < (2*m_fNumofChains); i++)
    {
        tree < CXYVector <float> >::iterator TempI;
        t1.push_back(TempI);
        r1.push_back(TempI);
    }
  //  top_r.push_back(t1);
  //  root_r.push_back(r1);
    top.push_back(t1);
    root.push_back(r1);
}

for (int i=0; i<m_samples; i++)
{
    for (int j=0; j<(2*m_fNumofChains); j++)
    {
       // top_r[i][j]= ptr_right[i][j]-> begin();
        top[i][j]= ptr[i][j]-> begin();

    }
}


int ind = 0;
CXYVector<float>  kV_centrpoint;
float t_alpha;
float t_alpha2;

//CXYVector<float> kV_endpoint;
//cout << "Size = " << m_pOctree.size() << endl;
for (int j=0; j < m_samples; j++)
{
    boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
    m_pOctree.push_back(new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator) );

    for ( int i = 0; i < (m_fNumofChains); i++ )
    {


        if (i==11)
        {
            kV_centrpoint = RndSetCentrPoint_nucleolus();
            t_alpha = getInitialAngle(kV_centrpoint, ind);
            t_alpha2 = getInitialAngle(kV_centrpoint, ind+1);
            cout<<"t_alpha is "<<t_alpha<<endl;

        }
        else
        {
            kV_centrpoint = RndSetCentrPoint(ind);
            t_alpha = getInitialAngle(kV_centrpoint, ind);
            t_alpha2 = getInitialAngle(kV_centrpoint, ind+1);
            cout<<"t_alpha is "<<t_alpha<<endl;
        }

        root[j][ind] = ptr[j][ind]-> insert (top[j][ind], kV_centrpoint);
        root[j][ind+1] = ptr[j][ind+1]-> insert (top[j][ind+1], kV_centrpoint);
        centromeres[j][ind][0] = kV_centrpoint[0];
        centromeres[j][ind][1] = kV_centrpoint[1];
        centromeres[j][ind][2] = kV_centrpoint[2];
        centromeres[j][ind+1][0] = kV_centrpoint[0];
        centromeres[j][ind+1][1] = kV_centrpoint[1];
        centromeres[j][ind+1][2] = kV_centrpoint[2];
        alpha[ind] = t_alpha;
        alpha[ind+1] = t_alpha2;




       // root_r[j][i] = ptr_right[j][i]-> insert (top_r[j][i], kV_startpoint);

    // boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
    // p_octree = new OcTree<float,int> (m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
    // m_pOctree.push_back(p_octree);
        m_pOctree[j] ->addPoint(kV_centrpoint[0],kV_centrpoint[1],kV_centrpoint[2],1,m_maxDepth);
//        CXYVector<float> cpoint = RndSetEndPoint(ind,kV_centrpoint);
//        endpoints[i][ind][0] = cpoint[0];
//        endpoints[i][ind][1] = cpoint[1];
//        endpoints[i][ind][2] = cpoint[2];
//      //  cout<<"are we here"<<endl;
//        CXYVector<float> epoint = RndSetEndPoint(ind+1,kV_centrpoint);
//        endpoints[i][ind+1][0] = epoint[0];
//        endpoints[i][ind+1][1] = epoint[1];
//        endpoints[i][ind+1][2] = epoint[2];

     //cout<<"adress is "<<m_pOctree[i]<<endl;

        ind = ind+2;

   // m_pOctree.push_back(p_octree);
// cout<<"address"<<*(ptr_right)<<endl;
    }

    ind = 0;

}

//////////////////////////////////////// Im here! //////////////////////////////////////////////////////////////////




cout<<"initialize chain has ended"<<endl;

}

//---------------------------------------------------------------------------------------------------------
float CXYEnsemble::getInitialAngle(CXYVector <float> cpoint, int j)
{

    float alp;

    vector <int> cl = GetChainLengths();

    float Lp = GetPersistenceLength();

    float L = (cl[j]-1)*Lp;
  //  cout << "L is "<<L<<endl;

    float d = sqrt(pow(cpoint[0],2)+pow(cpoint[1],2)+pow(cpoint[2],2));
  //  cout << "d is "<<d<<endl;


    float R = GetNucleusSphereDiameter()/2;
  //  cout << "R is "<<R<<endl;

    float nL = d + R;

    if (L > nL)
    {
        alp = 0;
    }
    else
    {
        alp = acos((pow(d,2)+pow(L,2)-pow(R,2))/(2*d*L)) * 180.0 / 3.14;
    }

    return alp;
}

//----------------------------------------------GAMZE 10/2-------------------------------------------
vector < vector < tree< CXYVector<float> >* >  > CXYEnsemble::GetTree()
{
    //cout<<"get tree"<<endl;
	return m_pTChain;
}


//------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------
CXYVector<float> CXYEnsemble::RndSetStartPoint(CXYVector<float>  kV_centrpoint)
{
	float fdiameter = GetNucleusSphereDiameter();
	float fradius = fdiameter/2;
	float nucleolus = 3000;
	MTRand mtrand;

    float theta;
    float phi;
    float ro;
    float offset;

    double PI = 3.14159265;
	// Check if the random point located in the neuclus sphere.
	// When one point satisfy the condition, we accept it.
	float fCoord[3];
	while (1) {
	    theta = mtrand.randExc(360);
	    phi = mtrand.randExc(360);
	    offset = mtrand.randExc(nucleolus);
	    ro = fradius - offset;
		fCoord[0] = ro * cos(theta * PI / 180.0) * sin(phi * PI / 180.0);
		fCoord[1] = ro * sin(theta * PI / 180.0) * sin(phi * PI / 180.0);
		fCoord[2] = ro * cos(theta * PI / 180.0);
		float dist = sqrt(pow((fCoord[0]-kV_centrpoint[0]),2)+pow((fCoord[1]-kV_centrpoint[1]),2)+pow((fCoord[2]-kV_centrpoint[2]),2));
	//	if (IsOnNE(fCoord) && fCoord[0]<=500 && dist<11000){ GG_28_1
	if (IsOnNE(fCoord) && dist<11000){
			break;
		}
	}
//	cout << x << ";" << y << ";" << z<<";" <<endl;
//	cout << x*x+y*y+z*z << endl;
//	cout << fradius2 << endl;
	CXYVector<float> kV (3,fCoord);
	return kV;
}
//----------------------------------------------------------------------------
CXYVector<float> CXYEnsemble::RndSetEndPoint(int j, CXYVector <float> cpoint)
{
	float fdiameter = GetNucleusSphereDiameter();
	float fradius = fdiameter/2;
	float nucleolus = 3000;
	int nofc= GetNumChains();
	CXYVector <int> centr = getCentromeres(nofc);
	vector <int> clengths = GetChainLengths();
    //cout<<"clength[j] is "<<clengths[j]<<endl;

	MTRand mtrand;


    float theta;
    float phi;
    float ro;
    float offset;
    float dist;

    double PI = 3.14159265;
	// Check if the random point located in the neuclus sphere.
	// When one point satisfy the condition, we accept it.
	float fCoord[3];
	while (1) {
	    theta = mtrand.randExc(360);
	    phi = mtrand.randExc(360);
	    offset = mtrand.randExc(nucleolus);
	    ro = fradius;
		fCoord[0] = ro * cos(theta * PI / 180.0) * sin(phi * PI / 180.0);
		fCoord[1] = ro * sin(theta * PI / 180.0) * sin(phi * PI / 180.0);
		fCoord[2] = ro * cos(theta * PI / 180.0);
		CXYVector<float> kV_temp (3,fCoord);
		dist = sqrt ( pow((fCoord[0]-cpoint[0]),2) + pow( (fCoord[1]-cpoint[1]),2)  + pow( (fCoord[2] - cpoint[2]),2));
	//	cout<<"dist is "<<dist<<endl;
		float thres = (clengths[j]*1500);
	//	cout<<"threshold is "<<thres<<endl;
	//	if (IsOnNE(fCoord) && fCoord[0]<=500 && dist<=thres){ GG_28_1
	//	    cout<<"hello"<<endl;
	if (IsOnNE(fCoord) && dist<=thres){
			break;
		}
	}
//	cout << x << ";" << y << ";" << z<<";" <<endl;
//	cout << x*x+y*y+z*z << endl;
//	cout << fradius2 << endl;

	CXYVector<float> kV (3,fCoord);
  //  m_pOctree[j] ->addPoint(kV[0],kV[1],kV[2],1,m_maxDepth);

	return kV;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
CXYVector<float> CXYEnsemble::RndSetCentrPoint(int j)
{
	float fdiameter = GetSBPSphereDiameter();
	float fradius = fdiameter/2;
	float nucleolus = 3000;
		int nofc= GetNumChains();
		float dist;

		CXYVector <int> centr = getCentromeres(nofc);

	MTRand mtrand;


    float x;
    float y;
    float z;
    float offset;

    double PI = 3.14159265;
	// Check if the random point located in the neuclus sphere.
	// When one point satisfy the condition, we accept it.
	float fCoord[3];
	while (1) {
	    x = mtrand.randExc(10000)-mtrand.randExc(10000);
	    y = mtrand.randExc(fradius)-mtrand.randExc(fradius);
	    z = mtrand.randExc(fradius)-mtrand.randExc(fradius);
	    //ro = fradius - offset;
		fCoord[0] = x;
		fCoord[1] = y;
		fCoord[2] = z;
		CXYVector<float> kV_temp (3,fCoord);
		//		dist = sqrt ( pow((fCoord[0]-epoint[0]),2) + pow( (fCoord[1]-epoint[1]),2)  + pow( (fCoord[2] - epoint[2]),2));

	//	cout<<x<<" "<<y<<" "<<z<<" "<<endl;
		if (IsInsideSBPSphere(fCoord) &&  fCoord[0]<=3000 && IsEqualityHolds(kV_temp,j)){
//	if (IsInsideSBPSphere(fCoord)){
			break;
		}
	}
	cout << x << ";" << y << ";" << z<<";" <<endl;
//	cout << x*x+y*y+z*z << endl;
//	cout << fradius2 << endl;

	CXYVector<float> kV (3,fCoord);
 //   m_pOctree[j] ->addPoint(kV[0],kV[1],kV[2],1,m_maxDepth);

	return kV;
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsEqualityHolds(CXYVector<float> point, int j)
{
    vector <int> cl = GetChainLengths();
    bool Flag;
    float R = GetNucleusSphereDiameter()/2;
    float d = sqrt(pow(point[0],2)+pow(point[1],2)+pow(point[2],2));
    float nR = d+(cl[j]-1)*1500;
    cout<< "R is "<<R<<endl;
    cout<<"nR is "<<nR<<endl;


    if (R > nR)
    {
        Flag = false;
    }
    else
    {
        Flag = true;
    }

    return Flag;
}

//------------------------------------------------------------------------------------

CXYVector<float> CXYEnsemble::RndSetCentrPoint_nucleolus(void)
{
	float fdiameter = GetSBPSphereDiameter();
	float fradius = fdiameter/2;
	float nucleolus = 3000;
		int nofc= GetNumChains();
		float dist;

		CXYVector <int> centr = getCentromeres(nofc);

	MTRand mtrand;


    float x;
    float y;
    float z;
    float offset;

    double PI = 3.14159265;
	// Check if the random point located in the neuclus sphere.
	// When one point satisfy the condition, we accept it.
	float fCoord[3];
	while (1) {
	    x = mtrand.randExc(3000)+1000;
	    y = mtrand.randExc(fradius)-mtrand.randExc(fradius);
	    z = mtrand.randExc(fradius)-mtrand.randExc(fradius);
	    //ro = fradius - offset;
		fCoord[0] = x;
		fCoord[1] = y;
		fCoord[2] = z;
		CXYVector<float> kV_temp (3,fCoord);
		//		dist = sqrt ( pow((fCoord[0]-epoint[0]),2) + pow( (fCoord[1]-epoint[1]),2)  + pow( (fCoord[2] - epoint[2]),2));

	//	cout<<x<<" "<<y<<" "<<z<<" "<<endl;
		if (fCoord[0]>=3000 && IsInsideSphere(kV_temp)){
	//if (IsInsideSBPSphere(fCoord)){
			break;
		}
	}
	cout << x << ";" << y << ";" << z<<";" <<endl;
//	cout << x*x+y*y+z*z << endl;
//	cout << fradius2 << endl;

	CXYVector<float> kV (3,fCoord);
 //   m_pOctree[j] ->addPoint(kV[0],kV[1],kV[2],1,m_maxDepth);

	return kV;
}
bool CXYEnsemble::IsOnNE(float* fCoord)
{
//    cout<<"girdik mi buraya"<<endl;
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	float offset = 500;
	float nucleolus = 3000;
	float fradius3 = (fradius-offset)*(fradius-offset);
	bool flag;
	if ((fCoord[0] <= nucleolus) && (fCoord[0]*fCoord[0] + fCoord[1]*fCoord[1] + fCoord[2]*fCoord[2]
			<= fradius2) && (fCoord[0]*fCoord[0] + fCoord[1]*fCoord[1] + fCoord[2]*fCoord[2]
			>= fradius3)) {
		flag = true;
	}else {
		flag = false;
	}
//	cout<<"ciktik mi burdan"<<endl;
	return flag;
}


//----------------------------------------------------------------------------
bool CXYEnsemble::IsInsideSphere(CXYVector<float> kV_point)
{

    //cout<< "girdik mi ki buraya"<<endl;
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	bool flag;

	if (kV_point.SquaredLength() < fradius2) {
	  //   cout<<"girdik mi buraya"<<endl;

		flag = true;

	}else {
		flag = false;
	}

    //cout<<"ciktik burdan"<<endl;
   // cout<<"flag is "<<flag<<endl;
	return flag;
	//cout<<"ciktik burdan"<<endl;


}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsInsideSphere(float* fCoord)
{
//    cout<<"girdik mi buraya"<<endl;
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	bool flag;
	if (fCoord[0]*fCoord[0] + fCoord[1]*fCoord[1] + fCoord[2]*fCoord[2]
			< fradius2) {
		flag = true;
	}else {
		flag = false;
	}
//	cout<<"ciktik mi burdan"<<endl;
	return flag;
}

//---------------------------------------GAMZE 10/11-------------------------------------
bool CXYEnsemble::IsInsideSBPSphere(float* fCoord)
{
//    cout<<"girdik mi buraya"<<endl;
	float fradius = GetSBPSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	float SBPCenter[3];
	SBPCenter[0] = -7000;
	SBPCenter[1] = 0;
	SBPCenter[2] = 0;
	float dist = sqrt (pow((fCoord[0]-SBPCenter[0]),2) + pow((fCoord[1]-SBPCenter[1]),2) + pow((fCoord[2]-SBPCenter[2]),2) );
	bool flag;
	if (dist<=fradius) {
		flag = true;
	}else {
		flag = false;
	}
//	cout<<"ciktik mi burdan"<<endl;
	return flag;
}
//-------------------------------------------------------------------------------------------------
bool CXYEnsemble::IsCollision(CXYVector<float> kV_point, int j)
{

   OcTree<float,int> * p_octree;

  float x_ = kV_point[0];
  float y_ = kV_point[1];
  float z_ = kV_point[2];
  //cout<<x_<<" "<<y_<<" "<<z_<<endl;

  float min_x_ = max(x_ - m_dr, -m_fNucleusSphereDiameter);
  float min_y_ = max(y_ - m_dr, -m_fNucleusSphereDiameter);
  float min_z_ = max(z_ - m_dr, -m_fNucleusSphereDiameter);
  Eigen::Matrix<float,4,1> minPosition_ (min_x_,min_y_,min_z_,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+m_dr,y_+m_dr,z_+m_dr,1);
   vector< OcTreeNode< float, int >* > nodes_;


//cout<<"m_pOctree is"<<m_pOctree[j]<<endl;
//boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
 //    p_octree = new OcTree<float,int> (m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
 //    m_pOctree.push_back(p_octree);

   m_pOctree[j]->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,m_maxDepth,true);

  // delete p_octree;

 //   cout<<"this is working"<<endl;


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
	//return false;
  }else {
    return false;
  }
}


//--------------------------------GAMZE 10/3------------------------------------------
bool CXYEnsemble::GrowmChain(int k, int f, int m_samples)
{




     bool Flag=true;
//-----------initialize m chain if it's already not initialized---------------------------------------------------------

CXYVector<int> centr(m_fNumofChains);
centr = getCentromeres(m_fNumofChains);
//--------------------------------------------------------------------------------------------



	MTRand mtrand;

	int iNumNodes = GetNumNodes();
	vector <int> ChainLenghts = GetChainLengths();

    CXYVector<int> armlenghts(m_fNumofChains);

	for (int i=0; i<m_fNumofChains; i++)
	{
	    armlenghts[i] = ChainLenghts[i] - centr[i];
	}



//----------------------------vector of matrices---------------------------------------------
	vector <vector <tree < CXYVector <float> >* > > ptr = GetTree();


//---------------------------------------------------------------------------------------------

if (f==0)
{
    InitializeChain(m_samples);
}




//-----------------------------initialize the weights----------------------------------


//---------------------------------------------------------------------------------------------
//  iNumNodes = 3;

//-------------------grow m chain-------------------------------------------------------
//-------------------each chain has k number of nodes------------------------
//two scenario: 1) if chain is already grown means f>0----------------------------------
//--------------------2) we just started to grow means f=0 ------------------------------------


    if (f>0)
        {
            int a=f-1;

   // int numnodes=(f*k) +k;
 //   int numnodes = centromeres[]
    cout<<"f is "<<f<<endl;

    int dir;
    int j=0;

    CXYVector< float> points(3);
    int counter=0;

     int chainid=0;

    int index[2*m_fNumofChains];

    for (int kk=0; kk<2*m_fNumofChains; kk++)
    {
        index[kk]=kk;
    }

    random_shuffle(index, index+2*m_fNumofChains);

    for (int i=0; i<m_samples; i++)
    {


     //   cout<<"this is happenning"<<endl;
        while(j<2*m_fNumofChains)
        {
            chainid = index[j];
            cout<<"chainid id is "<<chainid<<endl;
            cout<<"chainlength is "<<ChainLenghts[chainid]<<endl;
            int numnodes=2;
            points[0] = endpoints[i][chainid][0];
            points[1] = endpoints[i][chainid][1];
            points[3] = endpoints[i][chainid][2];

            while (numnodes<=ChainLenghts[chainid])

            {

               cout<<"numnodes is "<<numnodes<<endl;
               counter=counter+1;
               if(counter>500)
                    {
                        cout<<"counter is "<<counter<<endl;
                        exit(0);
                    }



         //       counter=0;
                if ((ChainLenghts[chainid]-numnodes)>0)
                {
                    dir = 0;
                }
                if (numnodes == ChainLenghts[chainid])
                {
                    dir = 1;
                }

                if (chainid==23 && (numnodes == 1 || numnodes == 2 || numnodes == 3 || numnodes == 4 ))
                {
                    dir = 2;
                }

                if (chainid==23 && numnodes > 4 && numnodes < ChainLenghts[chainid])
                {
                    dir = 3;
                }

                if (chainid==23 && numnodes == ChainLenghts[chainid])
                {
                    dir = 1;
                }

                if (chainid==22 && numnodes < 11)
                {
                    dir = 4;
                }

                if (chainid==22 && numnodes == 11)
                {
                    dir = 5;
                }

                if(GrowOneChain(i,chainid,numnodes,a,dir,points))
                {
                    //j=j+1;
                    numnodes=numnodes+1;
                    a=a+1;
                    counter = 0;
         //           <<"a is "<<a<<endl;
                }
                else
                {
                    numnodes=numnodes;
                    a=a;
                    counter= counter+1;

                }

          //      counter = 0;

//                else
//                {
//                    //j=j;
//                    counter=counter+1;
//                    if (counter >=5000)
//                    {
//                        indicator=1;
//                       // Flag = true;
//			//ms_logweight[i].first=0;
//                	    cout<<"counter is "<<counter<<endl;
//    	                pro_samp1[ind1]=i;
//        	            counter=0;
//                        j=j+1;
//                        ind1 = ind1+1;
//
//                    }
//                    cout<<"counter is "<<counter<<endl;
//                }
//                counter = 0;

            }
       //     else
        //    {
               j=j+1;
               a=0;
               cout<<"j is "<<j<<endl;
           }
           j=0;
        }
      //  j=0;
  //  }

//cout<<"is it out of loop"<<endl;


        }

 //cout<<"here I am"<<endl;
  return Flag;

}

bool CXYEnsemble::GrowOneChain(int m_samples, int j, int numnodes, int a, int dir, CXYVector<float> points)

{

        bool Flag=true;


        vector <vector< tree< CXYVector <float> >* > > ptr = GetTree();


        float sum = 0;


        //    CXYMatrix<float> points = GetPoints(ptr[m_samples][j],numnodes-1);
                  MTRand mtrand;


//	int iNumNodes = GetNumNodes();

     //   counter = 0;
        for ( int i = 1+ a; i < numnodes; i++ )
        {
  //          cout<<"good morning"<<endl;
            tree< CXYVector<float> >::pre_order_iterator node;
            node = ptr[m_samples][j] -> end();
            node --;




            CXYMatrix < float > kMSamplesPoints = GetSamplesOrg();
            kMSamplesPoints *= GetSegLength(i-1);



            vector<int> GoodPointInd;
            vector<int> NoCollisionPointInd;

       //     cout<<"node is "<<(*node)<<endl;

            GetGoodPoints( (*node), kMSamplesPoints, GoodPointInd, NoCollisionPointInd,m_samples, dir, i);

      //      cout<<" get good points HAPPENED"<<endl;
            int GoodPointSize = GoodPointInd.size();
            int NoCollisionPointSize = NoCollisionPointInd.size();
      //      cout<<"good point size is "<<GoodPointSize<<endl;
       //     cout<<"no collision point size is "<<NoCollisionPointSize<<endl;
            if (GoodPointSize == 0 || NoCollisionPointSize == 0) {
                Flag = false;
             //   counter = counter + 1;
                break;
            }


            CXYVector <float> beta(GoodPointSize);

            CXYVector <float> cent(3);

            int iRnd;
            int excellentpoints;
             CXYVector <int> excellentind(GoodPointSize);


                 CXYVector <int> indi(GoodPointSize);
                 for (int jj=0; jj<GoodPointSize; jj++)
                 {
                     indi[jj] = GoodPointInd[jj];
          //           cout<<"indi is "<<indi[jj]<<endl;
                 }
                 //indi = GoodPointInd;
                 CXYVector <float> pr(GoodPointSize);
                 pr = CalculateProbForGrowth((*node), kMSamplesPoints,GoodPointSize, indi, j, dir, i, points);


               //  excellentpoints = 0;
                 int aa = 0;
                 for (int ii=0; ii<GoodPointSize; ii++)
                 {

                     if ( pr[ii] == 1)
                     {
       //                   cout<<"here"<<endl;
        //                  cout<<"good point ind is "<<GoodPointInd[ii]<<endl;
                    //     excellentpoints++;
                         excellentind[aa] = GoodPointInd[ii];
        //                cout<<"excellent ind is "<<excellentind[aa]<<endl;
                         aa = aa + 1;
                     }
                 }

                  vector <FloatIntPair> p;
//
//
//
   for (int kk = 0; kk<GoodPointSize; kk++)
    {
       p.push_back(FloatIntPair(pr[kk],GoodPointInd[kk]));
    }

    sort(p.begin(), p.end(), FloatIntPairCompare());

                 excellentpoints = aa;

float tmp;


int cons = 0;



    tmp= p[GoodPointSize-1].first;
    cout <<"tmp is "<<tmp<<endl;
    if ( tmp == 1)
    {
        cons = 1;
    }



//cout<<"constraint is "<<constraint<<endl;

            if (excellentpoints == 0) {
                    iRnd = p[GoodPointSize-1].second;
           }
           else
           {
                int nRnd = mtrand.randInt(excellentpoints-1);
//
                iRnd = excellentind[nRnd];
           }



         cout<<"iRnd is "<<iRnd<<endl;

           cout<<"excellentpoints are "<<excellentpoints<<endl;
//



//             if (dir == 0 || dir == 1 || dir == 3 || dir == 4 || dir == 5)
//             {
//
//                 cout<<"alpha is "<<alpha[j]<<endl;
//
//                if (alpha[j] != 0)
//                {
//
//                    cent[0] = centromeres[0][j][0];
//                    cent[1] = centromeres[0][j][1];
//                    cent[2] = centromeres[0][j][2];
//
//                    beta = CalculateBeta((*node), kMSamplesPoints,GoodPointSize, GoodPointInd, j, dir, i, cent);
//
//
//
//                    excellentpoints = getExcellentPoints(beta, j, GoodPointSize);
//
//                    excellentind = getExcellentIndexes(beta, j, GoodPointSize, GoodPointInd, excellentpoints);
//
//                    cout<<"good points are "<<GoodPointSize<<endl;
//                    cout<<"excell points are "<<excellentpoints<<endl;
//

//                }
//                else
//                {
//                    for (int ss=0; ss<GoodPointSize; ss++)
//                    {
//                        beta[ss] = 0;
//                    }
//
//                    excellentpoints = getExcellentPoints(beta, j, GoodPointSize);
//
//                    excellentind = getExcellentIndexes(beta, j, GoodPointSize, GoodPointInd, excellentpoints);
//
//                    cout<<"good points are "<<GoodPointSize<<endl;
//                    cout<<"excell points are "<<excellentpoints<<endl;
//
//                }
//             }
//             if ( dir == 2 )
//             {
//                 excellentpoints = getExcellentPoints_nucleolus((*node), kMSamplesPoints,GoodPointSize, GoodPointInd);
//
//                 excellentind = getExcellentIndexes_nucleolus((*node), kMSamplesPoints,GoodPointSize, GoodPointInd, excellentpoints);
//
//                cout<<"good points are "<<GoodPointSize<<endl;
//                cout<<"excell points are "<<excellentpoints<<endl;
//          //      int nRnd = mtrand.randInt(excellentpoints-1);
//           //     iRnd = excellentind[nRnd];
//             }
//

//
//         //   cout<<"excell points are "<<excellentpoints<<endl;
//
//            CXYVector <float> prob_for_growth(GoodPointSize);              //calculate beta
//            prob_for_growth = CalculateProbForGrowth((*node), kMSamplesPoints,excellentpoints, excellentind, j, dir, i, points);
//      //       int iRnd;
////
//             vector <FloatIntPair> p;
//
//
//
//    for (int kk = 0; kk<excellentpoints; kk++)
//    {
//        p.push_back(FloatIntPair(0,1));
//    }
////int ss=0;
////int mm=0;
////
//    for (int kk = 0; kk<excellentpoints; kk++)
//   {
////        if  (prob_for_growth[k]==1)
////        {
////            mm=mm+1;
////        }
//            p[kk].first = prob_for_growth[kk];
//            p[kk].second = kk;
////            p[ss].second =k;
////            ss=ss+1;
////
//    }
//
//
////
////    if (mm == 0 && (dir == 0 || dir ==3 || dir ==4))
////    {
////         Flag = false;
////             //   counter = counter + 1;
////         break;
////    }
////
//    sort(p.begin(), p.end(), FloatIntPairCompare());
////
//
//
//    for (int kk=1; kk<excellentpoints; kk++)
//    {
//        p[kk].first=p[kk].first+p[kk-1].first;
//    }
//
//
//    int a_rnd=mtrand.randInt(100);


//        float a_rnd=  static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/p[ss-1].first));
//
////         int a_rnd=0;
////        while (a_rnd == 0)
////         {
////                a_rnd = static_cast <int> (rand()) / (static_cast <int> (RAND_MAX/100));
////         }
//       cout<<"a_rnd is "<<a_rnd<<endl;
//       // int j = rand() % m_sample;
//      //  cout<<"j is "<<j<<endl;
//      //  cout<<"prob[j] is "<<prob[j]<<endl;
//      float co=0;
////      for (int k=0; k<GoodPointSize; k++)
////      {
////          if (prob_for_growth[k]==1)
////          {
////              iRnd = k;
////              co = co+1;
////              break;
////          }
////      }
//      if (co ==0)
//      {
//        int an;
//      for (int kk=0; kk<excellentpoints; kk++)
//      {
//        if (p[kk].first > (p[excellentpoints-1].first/a_rnd)  )
//        {
//            iRnd = p[kk].second;
//            cout<<"iRnd is "<<iRnd<<endl;
//            cout<<"p[kk].first is "<<p[kk].first<<endl;
//            break;
//        }
//      }
////      }




          // cout<<"logweight is"<<m_LogWeight[j].first<<endl;
     //       m_LogWeight[j].first += log(NoCollisionPointSize);
      //      m_vLogWeight.push_back(m_LogWeight[j].first);
           // sum = sum + m_LogWeight[j].first;
           int NofC=GetNumChains();
           CXYVector<int> centr = getCentromeres(NofC);

            CXYVector<float> kV_Point(3);
            CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);

            kV_Point = kMSamplesPoints.GetRow(iRnd) + (*node);

            GetMiddlePoints((*node), kV_Point, kM_MiddlePoints);

//            cout<<"angle is "<<alpha[j]<<endl;


//            if (alpha[j] == 0)
//            {
//
//                centromeres[0][j][0] = kV_Point[0];
//                centromeres[0][j][1] = kV_Point[1];
//                centromeres[0][j][2] = kV_Point[2];
//                alpha[j] = getNewAngle(kV_Point, j, i);
//
//             }
//
//
//            float temp_weight = -log(excellentpoints) - log(prob_for_growth[an]);
//
//            cout<<"prob_for_growth is "<<prob_for_growth[an]<<endl;

            ms_logweight[m_samples].first += log(tmp);

            cout<<"logweight is "<<ms_logweight[m_samples].first<<endl;




            ptr[m_samples][j] -> append_child(node, kV_Point);



            for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
                m_pOctree[m_samples] -> addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
            }

            //cout<<"m_pOctree is "<<m_pOctree[j]<<endl;

        }
   //     ms_logweight[m_samples].first + = m_LogWeight[j].first


        return Flag;

}
//----------------------------------------------------------------------------------------------------------------------------------
float CXYEnsemble:: getNewAngle(CXYVector <float> kV_Point, int chainid, int ind)
{
    float nangle;

    vector <int> cl = GetChainLengths();

    float Lp = GetPersistenceLength();

    float L = (cl[chainid] - ind - 1) * Lp;

    float d = sqrt(pow(kV_Point[0],2)+pow(kV_Point[1],2)+pow(kV_Point[2],2));

    float R =  GetNucleusSphereDiameter()/2;

    float nL = d + R;

    if (L > nL)
    {
        nangle = 0;
    }
    else
    {
        nangle = acos ( (pow(d,2) + pow(L,2) - pow(R,2)) / (2*d*L) ) * 180 / 3.14;
    }

    cout<<"nangle is "<<nangle<<endl;


    return nangle;
}

//----------------------------------------------
int CXYEnsemble::getExcellentPoints(CXYVector <float> beta, int j, int GoodPointSize)
{

  //  cout<<"here we are"<<endl;

    int expoints = 0;

    for (int i=0; i<GoodPointSize; i++)
    {
        if (beta[i] >= alpha[j])
        {
            expoints = expoints+1;
        }
    }

    return expoints;


}
//-----------------------------------------------------------

//-------------------------------------------------------------

int CXYEnsemble::getExcellentPoints_nucleolus(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int GoodPointSize, vector<int> goodpointind)
{

  //  cout<<"here we are"<<endl;

    int expoints = 0;
    CXYVector <float> kV_Point(3);

    for (int i=0; i<GoodPointSize; i++)
    {
        kV_Point = kMSamplesPoints.GetRow(goodpointind[i]) + prvnode;
        if ( kV_Point[0] >=3000 || kV_Point[0] <=4000)
        {
            expoints=expoints+1;
        }
    }

    return expoints;


}
//-----------------------------------------------------------------------------------------------------------------------------------

CXYVector <int> CXYEnsemble::getExcellentIndexes_nucleolus(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int GoodPointSize, vector<int> goodpointind, int excellentpoints)
{

  //  cout<<"here we are"<<endl;

    CXYVector <int> ind(excellentpoints);
    CXYVector <float> kV_Point(3);
    int am = 0;

    for (int i=0; i<GoodPointSize; i++)
    {
        kV_Point = kMSamplesPoints.GetRow(goodpointind[i]) + prvnode;
        if ( kV_Point[0] >=3000 || kV_Point[0] <=4000)
        {
            ind[am] = goodpointind[i];
            am=am+1;
        }
    }

    return ind;


}
//-----------------------------------------------------------------------------------------------------------------------------------
CXYVector <int> CXYEnsemble::getExcellentIndexes(CXYVector <float> beta, int j, int GoodPointSize, vector <int> GoodPointInd, int excellentpoints)
{

  //  cout<<"here we are"<<endl;

    CXYVector <int> indexes(excellentpoints);
    int am = 0;

    for (int i=0; i<GoodPointSize; i++)
    {

        //cout<<"beta is "<<beta[i]<<endl;
       // cout<<"alpha is "<<alpha[j]<<endl;


        if (beta[i] >= alpha[j])
        {
          
            indexes[am] = GoodPointInd[i];
         //   cout<<indexes[a]<<endl;
            am=am+1;
        }
    }

   // cout<<"a is "<<a<<endl;

    return indexes;


}

//---------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
int CXYEnsemble::getMaxofArray(CXYVector<int> array_, int length_)
{
    int maxi = array_[0];

    for (int i=1; i<length_; i++)
    {
        if (array_[i]>maxi)
        {
            maxi = array_[i];
        }
    }

    return maxi;
}

//----------------------------------------------------------------------GAMZE 10/9----------------------------------------------

float CXYEnsemble::getMinofArray(CXYVector<float> array_, int length_)
{
    float maxi = array_[0];

    for (int i=1; i<length_; i++)
    {
        if (array_[i]<maxi)
        {
            maxi = array_[i];
        }
    }

    return maxi;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
void CXYEnsemble::Resampling(int k, int m_samples)
{

  //  std::srand ( unsigned ( std::time(0) ) );

    int iNumNodes = GetNumNodes();

    vector <float> probability;

    float sumofweights = 0;

    vector <int> list;

    int  index=0;

    int listsize=0;

    int rand1;
    CXYVector<int> centr(m_fNumofChains);
    centr = getCentromeres(m_fNumofChains);
    int maxcentr = getMaxofArray(centr,m_fNumofChains);

    CXYVector<int> rightarmlength(m_fNumofChains);

    for (int i=0; i<m_fNumofChains; i++)
    {
        rightarmlength[i] = m_fChainLengths[i] - centr[i];
    }

    int maxright = getMaxofArray(rightarmlength,m_fNumofChains);


 //   vector <vector <tree< CXYVector<float> >*> > ptr;

    
 //   int f=0;
   // int t = 0;
    // int t=k+(f*k);
    int m = iNumNodes;
    CXYVector<int> VectorofIndex;

    indicator=0;

    vector <vector< tree< CXYVector <float> >* > > ptr = GetTree();

      //      vector <vector< tree< CXYVector <float> >* > > ptr_left = GetTree_left();

int f=0;
    while (f==0 && GrowmChain(k,f,m_samples))
    {
        f=f+1;

    }

  //      cout<<"initialize chain is started"<<endl;
 //       InitializeChain_left(m_samples);
  //      InitializeChain_right(m_samples);
  //      cout<<"initialize chain is done"<<endl;
  //      f=f+1;
  //  }


int t=k+(f*k);
k = iNumNodes/2;

       while (t<=iNumNodes && GrowmChain(k,f,m_samples))
        {

            f=f+1;
            t=iNumNodes+1;
//            if (GrowmChain(k,f,m_samples))
//                {
//                         t=k+(f*k);
//
//                double ESS = CalculateEffSampleSize(m_samples);
//                cout<<"ESS is "<<ESS<<endl;
//                if (ESS <= 0.3*m_samples)
//                {
////
//
//                    cout<<"t is "<<t<<endl;
//                    VectorofIndex= getResampleIndex(m_samples);
//                    CreateNewPointTree(t,m_samples, VectorofIndex,1);
//                    CreateNewPointTree(t,m_samples, VectorofIndex,0);
//                    cout<<"we created new point tree"<<endl;
//                    CreateNewOctree(t,m_samples);
//                    cout<<"we created new octree"<<endl;
//                    ReconstructLogweight(m_samples, VectorofIndex);
//                    cout<<"we reconstructed logweights"<<endl;
//                   // if (t<=50)
//                   // {
//                        f=f+1;
//                    //     t=k+(f*k);
//                   // }
//                   // else
//                   // {
//                   //     f=f+0.5;
//                   // }
//              }
//              else
//             {
//                  f=f+1;
//
//            }
//
//
//
//            }

    }
}

CXYVector <float> CXYEnsemble::getResampleprob(int m_sample)
{
    CXYVector<float> prob(m_sample);
    cout<< "prob for weights is being created"<<endl;
    vector <FloatIntPair> weights;

    for (int i = 0; i<m_sample; i++)
    {
        weights.push_back(FloatIntPair(0,1));
    }
    for (int i = 0; i<m_sample; i++)
    {
        weights[i].first = ms_logweight[i].first;
        weights[i].second = ms_logweight[i].second;
    }

    sort(weights.begin(), weights.end(), FloatIntPairCompare());

    float max_prob = weights[m_sample-1].first;
    cout<<"max prob is "<<max_prob<<endl;

    for (int i=0; i<m_sample; i++)
    {
        prob[i] = exp(ms_logweight[i].first-max_prob);
    }

    return prob;
}
CXYVector <int> CXYEnsemble::getResampleIndex(int m_sample)
{
   CXYVector<int> new_vector(m_sample);
    CXYVector<float> prob = getResampleprob(m_sample);

    cout<<"get index is working"<<endl;
    int s=0;
    //int j=0;

    while(s<m_sample)
    {
        float a_rnd= (float)rand()/(float)RAND_MAX;
      //  cout<<"a_rnd is "<<a_rnd<<endl;
        int j = rand() % m_sample;
      //  cout<<"j is "<<j<<endl;
      //  cout<<"prob[j] is "<<prob[j]<<endl;
        if (prob[j] > a_rnd && prob[j]>0)
        {
            new_vector[s]=j;
            s=s+1;
        }
    }
    cout<<"s is "<<s<<endl;
    return new_vector;


}
double CXYEnsemble::CalculateEffSampleSize(int m)
{
   // cout<<"are we here yet "<<endl;
    vector <FloatIntPair> weights;
    for (int i = 0; i<m; i++)
    {
        weights.push_back(FloatIntPair(0,1));
    }

    for (int i = 0; i<m; i++)
    {
        weights[i].first = ms_logweight[i].first;
        weights[i].second =ms_logweight[i].second;
    }

    sort(weights.begin(), weights.end(), FloatIntPairCompare());

    float maxWeight= weights[m-1].first;
    cout<<"max weight is "<<maxWeight<<endl;
   //  cout<<"max weight is "<<m_LogWeight[3].first<<endl;

    CXYVector <float> relativeWeight(m);
    double sumweight = 0;
    CXYVector <float> q(m);

    for (int i = 0; i < m; i++)
    {
      //  cout<<"are we here yet "<<endl;
    //    cout<<"weight is "<<m_LogWeight[i].first<<endl;
     //   cout<<"weight is "<<maxWeight<<endl;
        relativeWeight[i] = exp(ms_logweight[i].first-maxWeight);
        // relativeWeight[i] = 0;


        sumweight = sumweight + relativeWeight[i];

    }



    for (int i = 0; i < m; i++)
    {
        q[i] = relativeWeight[i] / sumweight;
    }

    double ESS = 0;

    for (int i = 0; i < m; i++)
    {
        ESS = ESS + pow(q[i],2);
    }

    ESS = 1 / ESS;

    return ESS;
}

//void CXYEnsemble::CreateNewPointTree(int k, int m_samples, CXYVector<int> vectorofI, int dir)
//{
//    vector <vector <vector < tree<CXYVector<float> >* > > >  ptr = GetTree();
//    if (dir ==0)
//    {
//
//
//        vector < vector < CXYMatrix<float> > > newpoints;
//        for (int i = 0; i < m_samples; i++)
//        {
//
//            vector < CXYMatrix <float> >  v;
//            for (int j=0; j<m_fNumofChains; j++)
//            {
//
//                CXYMatrix <float>  MatTemp;
//                v.push_back(MatTemp);
//            // delete TreeTemp;
//            }
//            newpoints.push_back(v);
//        // v.clear();
//
//        }
//
//
//
//   // int iNumNodes= k;
//        vector <  tree < CXYMatrix <float> >::iterator>  top, root;
//        cout<<"starting the iteration in createnewpointtree for left"<<endl;
//        for (int s=0; s<m_samples; s++)
//        {
//
//            for ( int j = 0; j < m_fNumofChains; j++ )
//            {
//                cout<<"is this happening?"<<endl;
//   //             cout<<"k is "<<k<<endl;
//                newpoints[s][j] = GetPoints(ptr[0][s][j],k) ;
//          //  cout<<"get points is happening"<<endl;
//            }
//
//        }
//
//        cout<<"get points ended"<<endl;
//
//
//
//
//
//        for (int s=0; s<m_samples; s++)
//        {
//
//            vector <  tree < CXYVector <float> >::iterator>  top, root;
//     //   cout<<"index is "<<ms_logweight[m_samples-s-1].second<<endl;
//            CXYMatrix<float> newnewpoints;
//
//            for (int i=0; i<m_fNumofChains; i++)
//            {
//                newnewpoints = newpoints[vectorofI[s]][i];
//
//    //     if ( !ptr.at(m_LogWeight[j].second) -> empty() )
//    //        {
//                ptr[0][s][i] -> clear();
//     //       }
//                top.push_back(ptr[0][s][i] -> begin());
//                float fCoord[3] = {newnewpoints[0][0],newnewpoints[0][1],newnewpoints[0][2]};
//                CXYVector<float> kV_startpoint (3,fCoord);
//                root.push_back((ptr[0][s][i] -> insert (top.at(i), kV_startpoint)));
//
//
//                for ( int j = 1; j < k; j++ )
//                {
//
//                    tree< CXYVector<float> >::pre_order_iterator node;
//
//                    node = ptr[0][s][i] -> end();
//                    node --;
//                    CXYVector<float> kV_Point(3);
//
//            //cout<<"i is "<<i<<endl;
//
//                    kV_Point[0] = newnewpoints[j][0];
//                    kV_Point[1] = newnewpoints[j][1];
//                    kV_Point[2] = newnewpoints[j][2];
//
//
//                    ptr[0][s][i] -> append_child(node, kV_Point);
//           //  cout<<"is this happening"<<endl;
//           //cout<<"i is "<<i<<endl;
//         //  cout<<"iNumnodes is "<<k<<endl;
//                }
//            }
//
//      //  newnewpoints.Deallocate();
//    }
//    }
//    else
//    {
//        vector < vector < CXYMatrix<float> > > newpoints;
//        for (int i = 0; i < m_samples; i++)
//        {
//
//            vector < CXYMatrix <float> >  v;
//            for (int j=0; j<m_fNumofChains; j++)
//            {
//
//                CXYMatrix <float>  MatTemp;
//                v.push_back(MatTemp);
//            // delete TreeTemp;
//            }
//            newpoints.push_back(v);
//        // v.clear();
//
//	}
//
//
//
//   // int iNumNodes= k;
//        vector <  tree < CXYMatrix <float> >::iterator>  top, root;
//        cout<<"starting the iteration in createnewpointtree for right"<<endl;
//        for (int s=0; s<m_samples; s++)
//        {
//
//            for ( int j = 0; j < m_fNumofChains; j++ )
//            {
//                newpoints[s][j] = GetPoints(ptr[1][s][j],k) ;
//          //  cout<<"get points is happening"<<endl;
//            }
//
//        }
//
//        cout<<"get points ended"<<endl;
//
//
//
//
//
//        for (int s=0; s<m_samples; s++)
//        {
//
//            vector <  tree < CXYVector <float> >::iterator>  top, root;
//     //   cout<<"index is "<<ms_logweight[m_samples-s-1].second<<endl;
//            CXYMatrix<float> newnewpoints;
//
//            for (int i=0; i<m_fNumofChains; i++)
//            {
//                newnewpoints = newpoints[vectorofI[s]][i];
//
//    //     if ( !ptr.at(m_LogWeight[j].second) -> empty() )
//    //        {
//                ptr[1][s][i] -> clear();
//     //       }
//                top.push_back(ptr[1][s][i] -> begin());
//                float fCoord[3] = {newnewpoints[0][0],newnewpoints[0][1],newnewpoints[0][2]};
//                CXYVector<float> kV_startpoint (3,fCoord);
//                root.push_back((ptr[1][s][i] -> insert (top.at(i), kV_startpoint)));
//
//
//                for ( int j = 1; j < k; j++ )
//                {
//
//                    tree< CXYVector<float> >::pre_order_iterator node;
//
//                    node = ptr[1][s][i] -> end();
//                    node --;
//                    CXYVector<float> kV_Point(3);
//
//            //cout<<"i is "<<i<<endl;
//
//                    kV_Point[0] = newnewpoints[j][0];
//                    kV_Point[1] = newnewpoints[j][1];
//                    kV_Point[2] = newnewpoints[j][2];
//
//
//                    ptr[1][s][i] -> append_child(node, kV_Point);
//           //  cout<<"is this happening"<<endl;
//           //cout<<"i is "<<i<<endl;
//         //  cout<<"iNumnodes is "<<k<<endl;
//                }
//            }
//
//      //  newnewpoints.Deallocate();
//    }
//    }
//
//
//cout<<"create new points ended"<<endl;
//
//    //newpoints.clear();
//
//   // newpoints.clear();
//}
//
//void CXYEnsemble::CreateNewOctree(int k, int m_samples)
//
//{
//
//
//    vector <vector <vector <tree<CXYVector<float> >* > > > ptr = GetTree();
//
//    OcTree<float,int>* p_octree;
//
//      for ( int i = 0; i < m_pOctree.size(); i++ )
//        {
//
//            //if (m_pOctree[i] -> root_) {
//               // cout<<"is this hapening"<<endl;
//                delete m_pOctree[i];
//
//            //}
//        }
//
//
//    m_pOctree.clear();
//
//    //int iNumNodes= k;
//        for (int i = 0; i < m_samples; i++)
//	{
//        boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
//        m_pOctree.push_back(new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator) );
//
//	}
//
//
//    for( int j=0; j<m_samples; j++)
//    {
//
// //       a = GetPoints(ptr[j], iNumNodes, j);
//
//        for (int s=0; s<m_fNumofChains; s++)
//        {
//
//        CXYMatrix<float> new_points_right(k,3) ;
//        new_points_right = GetPoints(ptr[1][j][s], k);
//
//        CXYMatrix<float> new_points_left(k,3) ;
//        new_points_left = GetPoints(ptr[0][j][s], k);
//
//
//
//         //   cout<<"new_points are "<<new_points[5][1]<<endl;
//     //       boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
//    // p_octree = new OcTree<float,int> (m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
//    // m_pOctree.push_back(p_octree);
//
//     m_pOctree[j] ->addPoint(new_points_right[0][1],new_points_right[0][2],new_points_right[0][3],1,m_maxDepth);
//     m_pOctree[j] ->addPoint(new_points_left[0][1],new_points_left[0][2],new_points_left[0][3],1,m_maxDepth);
//
//        for (int i=1; i<k; i++)
//
//        {
//
//
//            CXYMatrix < float > kMSamplesPoints = GetSamplesOrg();
//            kMSamplesPoints *= GetSegLength(i-1);
//            float fCoord_pre_right[3];
//            float fCoord_right[3];
//            float fCoord_pre_left[3];
//            float fCoord_left[3];
//            CXYVector<float> kV_Point(3);
//            CXYVector<float> kV_Point_pre(3);
//            CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);
//
//            for (int index=0; index<3; index++)
//            {
//                fCoord_pre_right[index]=new_points_right[i-1][index];
//                fCoord_right[index] = new_points_right[i][index];
//               // cout<<"kv points are "<<fCoord[index]<<endl;
//               fCoord_pre_left[index]=new_points_right[i-1][index];
//                fCoord_left[index] = new_points_left[i][index];
//            }
//
//
//            kV_Point = CXYVector<float> (3,fCoord_right);
//            kV_Point_pre = CXYVector<float> (3,fCoord_pre_right);
//
//            GetMiddlePoints(kV_Point_pre, kV_Point, kM_MiddlePoints);
//
//            for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
//                m_pOctree[j] -> addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
//            }
//
//            kV_Point = CXYVector<float> (3,fCoord_left);
//            kV_Point_pre = CXYVector<float> (3,fCoord_pre_left);
//
//            GetMiddlePoints(kV_Point_pre, kV_Point, kM_MiddlePoints);
//
//            for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
//                m_pOctree[j] -> addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
//            }
//
//                    //    cout<<"is this happening"<<endl;
//              //      new_points.Deallocate();
//                          //  delete p_octree;
//        }
//
//        }
//
//        //cout<<"is this happening"<<endl;
//   //     cout<<"j is "<<j<<endl;
//       // newpoints.clear();
//        //delete[] newnewpoints;
//
//
//
//    }
//
//
//}
//
//CXYMatrix<float>  CXYEnsemble::GetPoints(tree< CXYVector<float> >*  ptr, int k)
//{
//
//    tree < CXYVector <float> > ::iterator pos;
//    CXYMatrix<float> oldpoints(k,3);
//    cout<<"k is "<<k<<endl;
//    cout<<"inside the get points"<<endl;
// //   vector < CXYMatrix<float> > a;
//
// //   for (int s=0; s<m; s++)
//  //  {
//    pos = ptr -> begin();
//    int i=0;
//    while(ptr->is_valid(pos))
//    {
//     //   cout<<"is this happening"<<endl;
//        oldpoints[i][0]=(*pos)[0];
//        oldpoints[i][1]=(*pos)[1];
//        oldpoints[i][2]=(*pos)[2];
//        cout<<"oldpoints are "<<oldpoints[i][2]<<endl;
//        ++pos;
//        i=i+1;
//        cout<<"i is "<<i<<endl;
//    }
//
//   return oldpoints;
//
// //   }
//
// //  cout<<"get points is ended"<<endl;
//
//}

void CXYEnsemble::ReconstructLogweight(int m, CXYVector<int> vectorofI)
{
   CXYVector<double> newweights(m);
    CXYVector<double> temp_weight(m);
    CXYVector<float> prob = getResampleprob(m);


    for (int i=0; i<m; i++)
    {
        newweights[i] = ms_logweight[i].first;
    }



    for (int i=0; i<m; i++)
    {
        temp_weight[i] = newweights[vectorofI[i]]-(log(prob[vectorofI[i]]));
      //  cout<<"log of prob is "<<log(2*prob[VectorofIndex[i]])<<endl;
    }

    for (int i=0; i<m; i++)
    {
	ms_logweight[i].first = temp_weight[i];
    }

   // newweights.Deallocate();
}
//---------------------here I am---------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void CXYEnsemble::WritemChain(char* fn, int m, int m_samples)
{
////http://stdcxx.apache.org/doc/stdlibug/34-2.html
//
	vector <vector < tree<CXYVector<float> >* > >  ptr = GetTree();

	vector <int> ChainLengths = GetChainLengths();

	int numnodes= GetNumNodes();

    float points[2*m_fNumofChains][numnodes][3];



// int count =0;


	char buffer[1024];

for (int j=0; j<m_samples; j++)
{

    ostream* fp;
        if (strcmp(fn,"")) {
        sprintf(buffer, "%s/%s_%d", GetOutPath(),fn,j);
        cout<<buffer<<endl;
		fp = new std::ofstream(buffer);
        }else {
            fp = &std::cout;
        }



       *fp << "# LogWeight= " << ms_logweight[j].first << endl;

int ind = 0;

	for (int i=0; i <  m_fNumofChains; i++)

	{
	    int count = 0;
	     tree < CXYVector <float> >::iterator pos;


	    pos=ptr[j][ind]->begin();
        //cout<<"chain lengths are "<<m_fChainLengths[i]<<endl;
    //    while (ptr[j][i]->is_valid(pos)) {
         while (count<ChainLengths[ind]){

             points[ind][count][0] = (*pos)[0];
             points[ind][count][1] = (*pos)[1];
             points[ind][count][2] = (*pos)[2];

 //          sprintf(buffer,"%3e\t%3e\t%3e",(*pos)[0],(*pos)[1],(*pos)[2]);
        //cout<<m_LogWeight[m_LogWeight[i].second].first<<endl;


//           *fp << buffer <<endl;
		//cout<<buffer<<endl;
            ++pos;
            count=count+1;
	//	cout<<"pos is "<<count<<endl;
        }
        count = 0;
        pos=ptr[j][ind+1]->begin();
        while (count<ChainLengths[ind+1]){

            points[ind+1][count][0] = (*pos)[0];
             points[ind+1][count][1] = (*pos)[1];
             points[ind+1][count][2] = (*pos)[2];

   //        sprintf(buffer,"%3e\t%3e\t%3e",(*pos)[0],(*pos)[1],(*pos)[2]);
        //cout<<m_LogWeight[m_LogWeight[i].second].first<<endl;


    //        *fp << buffer <<endl;
		//cout<<buffer<<endl;
            ++pos;
            count=count+1;
	//	cout<<"pos is "<<count<<endl;
        }
        ind = ind+2;

	}

	ind = 0;
	for (int k=0; k<m_fNumofChains; k++)
	{
	    	    *fp<< "#  Chr "<<(k+1)<<endl;
	    	    for (int t=0; t<ChainLengths[ind]; t++)
	    	    {
	    	        sprintf(buffer,"%3e\t%3e\t%3e",points[ind][m_fChainLengths[ind]-t-1][0],points[ind][m_fChainLengths[ind]-t-1][1],points[ind][m_fChainLengths[ind]-t-1][2]);
	    	        *fp << buffer <<endl;
	    	    }
	    	    for (int t=1; t<ChainLengths[ind+1]; t++)
	    	    {
	    	        sprintf(buffer,"%3e\t%3e\t%3e",points[ind+1][t][0],points[ind+1][t][1],points[ind+1][t][2]);
	    	        *fp << buffer <<endl;
	    	    }

	    	    ind = ind+2;
	}




	if (fp != &std::cout)
		delete fp;

	}



    }
//----------------------------------------------------------------------------

void CXYEnsemble::SetOutPath(char* cPathName)
{
	CXYFile::MakeDirectory(cPathName,0755);
  cout << "make directory " << cPathName << endl;
	strcpy(m_cOutPath,cPathName);
}

//----------------------------------------------------------------------------
char* CXYEnsemble::GetOutPath(void)
{
	return m_cOutPath;
}
//----------------------------------------------------------------------------
// void CXYEnsemble::SetSegLengths(char* cStartEndFile){
// 	CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);
// 	for (int i= 0; i<kMStartEnd.GetRows(); i++) {
// 		int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
// 		m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
// 	}
// }
//----------------------------------------------------------------------------
void CXYEnsemble::SetSegLengths(char* cStartEndFile,const char* cMethod){
	CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);

	if (strcmp(cMethod, "AVG") == 0) {
		// average different node length
		float len = 0L;
		for (int i= 0; i<kMStartEnd.GetRows(); i++) {
			len = len + ( kMStartEnd[i][1] - kMStartEnd[i][0] + 1 ) * GetPackingDensity();
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
		cerr	<< "Please type Method" << endl;
		exit(-1);
	}
}
//----------------------------------------------------------------------------

vector<float> & CXYEnsemble::GetSegLengths(void){
	return m_vfSegLength;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetSegLength(int ind){
	return m_vfSegLength[ind];
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetSegLengths(void){
	for (int i= 0; i<GetNumNodes(); i++) {
		m_vfSegLength.push_back(GetPersistenceLength());
	}
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetContIndex(char* cContIndFile){
//	CXYVector<int> vtmp = CXYFile::ReadVectorInt(cContIndFile);
//	for (int i =0; i<vtmp.GetSize(); i++) {
//		// index start from 0
//		m_MContInd.push_back(vtmp[i]-1);
//	}
  m_MContInd = CXYFile::ReadMatrix(cContIndFile);
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetContIndex(void){
  m_MContInd.SetSize(m_iNumNodes,3);
  for (int i=0; i<m_iNumNodes; i++) {
    m_MContInd[i][0] = i+1;
    m_MContInd[i][1] = i+1;
    m_MContInd[i][2] = i+1;
  }
}
//----------------------------------------------------------------------------
void CXYEnsemble::GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints){

  float dx = EndPoint[0] - StartPoint[0];
  float dy = EndPoint[1] - StartPoint[1];
  float dz = EndPoint[2] - StartPoint[2];

  for (int i =0; i< m_iNumMiddleEndPoints; i++) {
    MiddleEndPoints[i][0] = StartPoint[0] + ((i+1)/(float)m_iNumMiddleEndPoints) * dx;
    MiddleEndPoints[i][1] = StartPoint[1] + ((i+1)/(float)m_iNumMiddleEndPoints) * dy;
    MiddleEndPoints[i][2] = StartPoint[2] + ((i+1)/(float)m_iNumMiddleEndPoints) * dz;
  }

}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsSatisfyCondition(CXYMatrix<float> & MiddleEndPoints, int j){
  bool flag = true;
  for (int i=0; i< m_iNumMiddleEndPoints; i++) {
    


    if (! IsInsideSphere(MiddleEndPoints.GetRow(i))  ) {
        //cout << "i is "<<i<<endl;

      flag = false;
      break;
    }

    

    if (IsCollision(MiddleEndPoints.GetRow(i),j)){
//        cout<<"this has happenned"<<endl;
 //     cout << MiddleEndPoints[m_iNumMiddleEndPoints-1][0] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][1] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][2] << endl;
      flag = false;
      break;
    }
  }
  //cout<<"flag is"<<flag<<endl;
  return flag;
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsSatisfyCentr(CXYVector<float> & kV_Point, int j, int dir, int i){
  bool flag = true;

  int numofchains = GetNumChains();
  vector <int> chainlengths = GetChainLengths();
  CXYVector <int> centromeres = getCentromeres(numofchains);

  CXYVector <int> SBPcenter = (-7000,0,0);
  float SBPradius = 3000;

if(j!=11)
{
  if (dir == 0)
  {
     if (kV_Point[0]>3000)
        {
            flag = false;
          //  break;
        }
  }

  if (dir == 1)
  {
      if (kV_Point[0]>3000)
      {
            flag = false;
          //  break;
      }
  }
}
else
{
     if (dir == 0)
    {
     if (i<=28 && i>=33 && kV_Point[0]>3000)
        {
            flag = false;
          //  break;
        }
    }
        if (dir == 1)
    {
     if (i<=28 && i>=33 && kV_Point[0]>3000)
        {
            flag = false;
          //  break;
        }
  }


}


  return flag;
}
//-----------------------------------------
bool CXYEnsemble::IsCloseSBP(CXYVector<float> & kV_Point, float SBPradius, CXYVector<int> & SBPcenter, int centr, int i){
  bool flag = false;

  float Lp = GetPersistenceLength();

  float sq_dist = sqrt(pow((kV_Point[0]-SBPcenter[0]),2) + pow((kV_Point[1]-SBPcenter[1]),2) + pow((kV_Point[2]-SBPcenter[2]),2)) - SBPradius;

  float compared_dist =  ((0.6*Lp*(centr - i)) );

  if (sq_dist < compared_dist)
  {
      flag = true;
    //  break;
  }

  return flag;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
bool CXYEnsemble::IsCloseNE(CXYVector<float> & kV_Point, int chainlength, int i){
  bool flag = false;

  float Lp = GetPersistenceLength();

  float radius = GetNucleusSphereDiameter()/2 - 500;

  float sq_dist = sqrt(pow((kV_Point[0]-0),2) + pow((kV_Point[1]-0),2) + pow((kV_Point[2]-0),2));

  float compared_dist =   radius - ((0.6*Lp*(chainlength-i)));

  if (sq_dist > compared_dist)
  {
      flag = true;
    //  break;
  }

  return flag;
}
//----------------------------------------------------------------
bool CXYEnsemble::IsCollision(CXYMatrix<float> & MiddleEndPoints, int j){
  bool flag = false;

  for (int i=0; i< m_iNumMiddleEndPoints; i++) {
  //    cout<<"this is is "<<i<<endl;
    if ( IsCollision(MiddleEndPoints.GetRow(i),j)){
          //  cout << MiddleEndPoints[m_iNumMiddleEndPoints-1][0] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][1] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][2] << endl;
      flag = true;
      break;
    }
  }
  return flag;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------

CXYVector <float> CXYEnsemble::CalculateBeta(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int GoodPointSize, vector<int> goodpointind, int chainid, int dir, int nodeid, CXYVector <float> cent)
{

    CXYVector<float> kV_Point(3);

    CXYVector <float> bt(GoodPointSize);
    //vector <int> cl = GetChainLengths();


    float R = GetNucleusSphereDiameter()/2;
    float d = sqrt(pow(cent[0],2)+pow(cent[1],2)+pow(cent[2],2));

    int iSampleSize = GoodPointSize;


    for (int i = 0; i < iSampleSize; i++) {
       //     cout<<"is it here"<<endl;
            kV_Point = kMSamplesPoints.GetRow(goodpointind[i]) + prvnode;
            float d_prime = sqrt(pow(kV_Point[0],2)+pow(kV_Point[1],2)+pow(kV_Point[2],2));
            float Lp = sqrt(pow((cent[0]-kV_Point[0]),2)+pow((cent[1]-kV_Point[1]),2)+pow((cent[2]-kV_Point[2]),2));

         //   cout <<"d_prime is "<<d_prime<<endl;
          //  cout<<"Lp is "<<Lp<<endl;
          //  cout<<"d is "<<d<<endl;
            bt[i] = acos((pow(d,2)+pow(Lp,2)-pow(d_prime,2))/(2*Lp*d)) * 180 / 3.14;
         //   cout <<"bt i is "<<bt[i]<<endl;
    }

    return bt;

}
 CXYVector <float> CXYEnsemble::CalculateProbForGrowth(CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints, int excelpointsize, CXYVector <int> excelpointind, int m_samples, int dir, int ind, CXYVector<float> points)
 {

    //           cout<<"excell point size is "<<excelpointsize<<endl;
          CXYVector<float> prob_for_growth(excelpointsize);
          CXYVector<float> dist(excelpointsize);
          CXYVector<float> h(excelpointsize);
          CXYVector<float> kV_Point(3);
          int iSampleSize = excelpointsize;
           int numofchr = GetNumChains();
        CXYVector<int> centr = getCentromeres(numofchr);
        vector<int> cl = GetChainLengths();
        float SBPradius = GetSBPSphereDiameter()/2;
        float SBPcenter[3];
        SBPcenter[0] = -7000;
        SBPcenter[1] = 0;
        SBPcenter[2] = 0;
        float Ncenter[3];
        Ncenter[0] = 0;
        Ncenter[1] = 0;
        Ncenter[2] = 0;
        float NEradius = GetNucleusSphereDiameter()/2 - 500;
        int sep[iSampleSize];
        float sum_dist = 0;

        for (int i = 0; i < iSampleSize; i++) {
            h[i]=0;
    }

        float h_temp;


        for (int i = 0; i < iSampleSize; i++) {
            kV_Point = kMSamplesPoints.GetRow(excelpointind[i]) + prvnode;
   //        cout<<"goodpointind is "<<excelpointind[i]<<endl;
    //       cout<<"i is "<<i<<endl;
    //       cout<<"dir is "<<dir<<endl;

            if (dir == 0)
            {
                sep[i] = cl[m_samples] ;
         //       dist[i] = sqrt(pow((kV_Point[0]-points[0]),2)+ pow((kV_Point[1]-points[1]),2) + pow((kV_Point[2]-points[2]),2));
                dist[i] = NEradius -  sqrt(pow((kV_Point[0]-Ncenter[0]),2)+ pow((kV_Point[1]-Ncenter[1]),2) + pow((kV_Point[2]-Ncenter[2]),2));

                sum_dist = sum_dist+dist[i];
            }
            else if (dir ==1)
            {
                dist[i] = NEradius -  sqrt(pow((kV_Point[0]-Ncenter[0]),2)+ pow((kV_Point[1]-Ncenter[1]),2) + pow((kV_Point[2]-Ncenter[2]),2));
                if (dist[i]<0)
                {
                    dist[i] = 0;
                }
                sum_dist = sum_dist+dist[i];
            }
            else if (dir ==2)
            {
                dist[i] = kV_Point[0]-3000;
                if (dist[i]>=0)
                {
                   // dist[i]=kV_Point[0]-500; //NEradius -  sqrt(pow((kV_Point[0]-Ncenter[0]),2)+ pow((kV_Point[1]-Ncenter[1]),2) + pow((kV_Point[2]-Ncenter[2]),2));
                   dist[i]=0;
                }

                sum_dist = sum_dist+dist[i];
            }
            else if (dir ==3)
            {
                dist[i] = NEradius -  sqrt(pow((kV_Point[0]-Ncenter[0]),2)+ pow((kV_Point[1]-Ncenter[1]),2) + pow((kV_Point[2]-Ncenter[2]),2));
                sep[i] = cl[m_samples] - ind;

            }
            else if (dir == 4)
            {
                dist[i] = SBPradius -  sqrt(pow((kV_Point[0]-SBPcenter[0]),2)+ pow((kV_Point[1]-SBPcenter[1]),2) + pow((kV_Point[2]-SBPcenter[2]),2));
                sep[i] = centr[m_samples] - ind;
            }
            else if (dir == 5)
            {
                 dist[i] = SBPradius -  sqrt(pow((kV_Point[0]-SBPcenter[0]),2)+ pow((kV_Point[1]-SBPcenter[1]),2) + pow((kV_Point[2]-SBPcenter[2]),2));
                 if (dist[i]<=0)
                 {
                     dist[i]=0;
                 }

            }
        }

        for (int i = 0; i < iSampleSize; i++) {
            if (dir == 0 || dir ==3 || dir ==4)
            {
                if (abs(dist[i]) <= pow((sep[i] * 1500),1))
                {
                    h[i] = 0;
                }
                else
                {
                    h[i] = abs(dist[i])/2;
                }
            }
            if (dir == 1 || dir == 2 || dir == 5)
            {
                h[i] = abs(dist[i]);
               //h[i]=0;
            }

        }


                        //                if (ind == centr[m_samples])
//                {
//                    if(IsInsideSBPSphere(kV_Point))
//                    {
//                        h_temp=0;
//                        h[i]=h[i]+h_temp;
//                    }
//                    else
//                    {
//                        h_temp=distance_score_SBP(kV_Point);
//                        h[i]=h[i]+h_temp;
//                    }
//                }
//                else
//                {
//                    h_temp = satisfy_loop_condition(kV_Point,points);
//                    h[i]=h[i]+temp;
//
//                }
//
//
//
//          //      sep[i] = centr[m_samples] - ind - 1;
//         //       cout<<"dist is "<<dist[i]<<endl;
//            }
//
//            if (dir == 1 && m_samples!=11)
//            {
//                sep[i] = cl[m_samples] - ind ;
//                  if(ind == cl[m_samples])
//                  {
//                      if(IsOnNE(kV_Point))
//                      {
//                          h_temp=0;
//                          h[i]=h[i]+h_temp;
//                      }
//                      else
//                      {
//                          h_temp=distance_score_NE(kV_Point);
//                          h[i]=h[i]+h_temp;
//                      }
//                  }
//     //             else if(ind==29 || ind ==30 || ind == 31 || ind==32)
//                  else
//                  {
//                        h_temp = satisfy_loop_condition(kV_Point,points);
//                        h[i]=h[i]+temp;
//                  }
//          //      sum_dist = sum_dist + dist[i];
//
//               // cout<<"sep is "<<sep[i]<<endl;
//            }
//
//
//           if (dir == 1 && m_samples==11)
//            {
//                sep[i] = cl[m_samples] - ind ;
//                  if(ind == cl[m_samples])
//                  {
//                      if(IsOnNE(kV_Point))
//                      {
//                          h_temp=0;
//                          h[i]=h[i]+h_temp;
//                      }
//                      else
//                      {
//                          h_temp=distance_score_NE(kV_Point);
//                          h[i]=h[i]+h_temp;
//                      }
//                  }
//                  else if(ind<29)
//                  {
//                     h_temp = satisfy_loop_condition(kV_Point,rDNApoints);
//                     h[i]=h[i]+temp;
//                  }
//                 else if(ind==29 || ind ==30 || ind == 31 || ind==32)
//                 {
//                     if (kV_Point[0]>=500)
//                     {
//                         h_temp=0;
//                          h[i]=h[i]+h_temp;
//                     }
//                     else
//                     {
//                         h_temp=distance_score_NO(kV_Point);
//                         h[i]=h[i]+h_temp;
//                     }
//                 }
//                  else
//                  {
//                        h_temp = satisfy_loop_condition(kV_Point,points);
//                        h[i]=h[i]+temp;;
//
//                  }
//          //      sum_dist = sum_dist + dist[i];
//
//
//               // cout<<"sep is "<<sep[i]<<endl;
//            }
//        }

    //    cout<<"sum of the dist is "<<sum_dist<<endl;

        float distmin = getMinofArray(dist,iSampleSize);

        float sum_prob = 0;
        float distave;


          //      cout<<"h[i] is "<<h[i]<<endl;
  for (int i = 0; i < iSampleSize; i++) {

      prob_for_growth[i] = exp(-h[i]/10000);

        //    cout<<"prob is "<<prob_for_growth[i]<<endl;

        }

 //       cout<<"sum  prob is "<<sum_prob<<endl;

    //  cout<<"is there a problem"<<endl;

    //    cout<<"prob is "<<prob_for_growth[5]<<endl;

  return prob_for_growth;


 }

   float CXYEnsemble::CalculateTargetDistribution(CXYVector<float> kV_Point, int dir, CXYVector<float> points)
   {
       float targetdist;
       float temp;
       float Ncenter[3];
       Ncenter[0] = 0;
        Ncenter[1] = 0;
        Ncenter[2] = 0;
        float NEradius = GetNucleusSphereDiameter()/2 - 500;
         float SBPradius = GetSBPSphereDiameter()/2;
        float SBPcenter[3];
        SBPcenter[0] = -7000;
        SBPcenter[1] = 0;
        SBPcenter[2] = 0;

       if(dir == 0 || dir == 3 || dir ==4)
       {
           targetdist =exp(-0);
           cout<<targetdist<<endl;
       }
       if (dir ==1)
       {

           temp =  NEradius -  sqrt(pow((kV_Point[0]-Ncenter[0]),2)+ pow((kV_Point[1]-Ncenter[1]),2) + pow((kV_Point[2]-Ncenter[2]),2));
           if (temp<=0)
           {
               targetdist=exp(0);
           }
           else
           {
                targetdist = exp(-temp/150);
           }
                      cout<<"temp is "<<temp<<endl;

       }
       if (dir == 2)
       {
            temp = (kV_Point[0]-3000);
            if (temp>0)
           {
               targetdist=exp(0);
           }
           else
           {
                targetdist = exp(-abs(temp)/150);
           }
                      cout<<"temp is "<<temp<<endl;

       }
       if (dir ==5)
       {
           temp = SBPradius -  sqrt(pow((kV_Point[0]-SBPcenter[0]),2)+ pow((kV_Point[1]-SBPcenter[1]),2) + pow((kV_Point[2]-SBPcenter[2]),2));
            if (temp<=0)
           {
               targetdist=exp(0);
           }
           else
           {
                targetdist = exp(-abs(temp)/150);
           }
                      cout<<"temp is "<<temp<<endl;
       }



       return targetdist;

   }

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
void CXYEnsemble::GetGoodPoints( CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints,vector<int>& GoodPointInd, vector<int>& NoCollisionPointInd, int j, int dir, int i){

  cout << prvnode[0] << " " << prvnode[1] << " " << prvnode[2] << endl;

  CXYVector<float> kV_Point(3);
  CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);

  int iSampleSize = kMSamplesPoints.GetRows();



//  cout << "===================="<< iSampleSize << endl;
  for (int i = 0; i < iSampleSize; i++) {
    kV_Point = kMSamplesPoints.GetRow(i) + prvnode;

   // cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;

    GetMiddlePoints(prvnode, kV_Point, kM_MiddlePoints);

 //   cout << "i = " << i <<endl;
 //   cout<<"km middle points is "<<kM_MiddlePoints[1][2]<<endl;
 //   cout<<"j is "<<j<<endl;
 if (dir == 2)
 {
     if (IsSatisfyCondition(kM_MiddlePoints, j) &&  (kV_Point[0]>=3000)){
    //    cout << "this part has happenned" <<endl;
      GoodPointInd.push_back(i);
    }
    // No collision
    if (! IsCollision(kM_MiddlePoints, j) &&  (kV_Point[0]>=3000)){
      NoCollisionPointInd.push_back(i);
    }
 }
 else
 {
    if (IsSatisfyCondition(kM_MiddlePoints, j) && (kV_Point[0]<=3000)){
    //    cout << "this part has happenned" <<endl;
    // if (IsSatisfyCondition(kM_MiddlePoints, j)){
      GoodPointInd.push_back(i);
    }
    // No collision
    if (! IsCollision(kM_MiddlePoints, j) && kV_Point[0]<=3000){
 //  if (! IsCollision(kM_MiddlePoints, j)){
      NoCollisionPointInd.push_back(i);
    }
 }
  }
//  cout << "===================="<< iSampleSize << endl;
    //cout<<"is this happening"<<endl;


}

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
