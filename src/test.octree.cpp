/* Test octree 
*/

#include <iostream>
#include "tree.hh"
#include "XYSO3Sequence.h"
#include "XYVector.h"
#include "octree.h"
#include <vector>

using namespace spatialaggregate;
int main (int argc, char *argv[]){
  Eigen::Matrix< float, 4, 1 > center(0.0f,0.0f,0.0f,1);
//  Eigen::Matrix< float, 4, 1 > center(0.1f,0.2f,0.3f,1);
  float minimumVolumeSize = 1; 
  float maxDistance = pow(2.0,8); 
  boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  = boost::make_shared< OcTreeNodeAllocator< float , int > >();
  OcTree<float,int> octree_(center,minimumVolumeSize,maxDistance,allocator);
  int maxDepth = ceil(octree_.depthForVolumeSize(1));
  cout << "maxDepth = " <<  maxDepth << endl;
  

  // Get position
//  OcTreeKey<float,int> octreekey_;
//  octreekey_ = octree_.getKey(2.2, 2.3, 12.4 );
//  Eigen::Matrix<float,4,1> pos_ = octreekey_.getPosition(&octree_);
//  cout << pos_[0] << "," << pos_[1] << "," << pos_[2]<< endl;
  
  // Add point
  OcTreeNode<float,int>* poctreenode_;

  poctreenode_ = octree_.addPoint(0, 0, 0, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 0, 0, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 1, 0, 1, maxDepth);
  poctreenode_ = octree_.addPoint(0, 1, 0, 1, maxDepth);
 
  poctreenode_ = octree_.addPoint(0, 0, 1, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 0, 1, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 1, 1, 1, maxDepth);
  poctreenode_ = octree_.addPoint(1, 1, 1.01, 1, maxDepth);
  poctreenode_ = octree_.addPoint(0, 1, 1, 1, maxDepth);

//  poctreenode_ = octree_.addPoint(2.5, 2.5,2.5, 1, maxDepth);
  
  
  
  
//  OcTreeNode<float,int>* prepresentativetarget_;
//  // algorithm O(log(n))
//  OcTreeNode<float,int>* prepresentative_ = octree_.findRepresentative(currentpos_,maxDepth);
//  // trace back to grand parent to search if the children are in the range.
//  if (prepresentative_) {
//    prepresentativetarget_ = prepresentative_;
//    OcTreeNode<float, int>* prepresentativeparent_ = prepresentative_->parent_;
//    if (prepresentativeparent_) {
//      prepresentativetarget_ = prepresentativeparent_;
//      OcTreeNode<float, int>* prepresentativegrandparent_ = prepresentativeparent_->parent_;
//      if (prepresentativegrandparent_) {
//        prepresentativetarget_ = prepresentativegrandparent_;
//      }
//    }
//  }else {
//    return false;
//  }
//  
//  float x_ = 0.5f;
//  float y_ = 0.5f;
//  float z_ = 0.5f;
//  float x_ = 2.0f;
//  float y_ = 2.5f;
//  float z_ = 2.5f;
  float x_ = 1.0f;
  float y_ = 1.0f;
  float z_ = 1.05f;
  float dr = pow(2.0,8);
  
  Eigen::Matrix<float,4,1> currentpos_(x_,y_,z_,1);
  cout << "Test position" << endl;
  cout << currentpos_[0] <<"," << currentpos_[1] << "," << currentpos_[2] << endl;
  cout << "dr = " << dr << endl;
  
  vector< OcTreeNode< float, int >* > nodes_;
  Eigen::Matrix<float,4,1> minPosition_ (x_-dr,y_-dr,z_-dr,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+dr,y_+dr,z_+dr,1);
//  prepresentativetarget_->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,maxDepth,true);
  //-------------------------------
  OcTreeKey<float,int> octreekey_;
  Eigen::Matrix<float,4,1> pos_;
  octreekey_ = octree_.getKey(minPosition_ );
  pos_ = octreekey_.getPosition(&octree_);
  cout << "minPosition_ " << pos_[0] << "," << pos_[1] << "," << pos_[2]<< endl;
  octreekey_ = octree_.getKey(maxPosition_ );
  pos_ = octreekey_.getPosition(&octree_);
  cout << "maxPosition_ " << pos_[0] << "," << pos_[1] << "," << pos_[2]<< endl;
  //-------------------------------
  
  octree_.getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,maxDepth,true);
  cout << "collision numbers = " << nodes_.size() << endl;
  
  if (nodes_.size() > 0) {
    cout << "-------this is collision-------" << endl;
  }else {
    cout << "-------no collision, we need add point.------" << endl;
  }

  for (unsigned int i=0; i<nodes_.size(); i++) {
    Eigen::Matrix<float,4,1> neighborpos_ = (nodes_[i])->getPosition();
    cout << nodes_[i] << "\t";
    cout << neighborpos_[0] << "," << neighborpos_[1] << "," << neighborpos_[2]  << endl;
  }




  //---------------debug---------------
    
  Eigen::Matrix< float, 4, 1 > m_center;
  m_center = Eigen::Matrix< float, 4, 1 >(0.0f,0.0f,0.0f,1);
  m_center = Eigen::Matrix< float, 4, 1 >(236.693, 995.11, 482.206, 1);
  
  float m_fCollisionLength = 100.0f;
  float m_fNucleusSphereDiameter = 6000.0f;
  float m_minimumVolumeSize = m_fCollisionLength; 
  float m_dr = m_minimumVolumeSize;
  allocator = boost::make_shared< OcTreeNodeAllocator< float , int > >();
  OcTree<float,int>* m_pOctree;
  m_pOctree = new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
  int m_maxDepth = ceil(m_pOctree->depthForVolumeSize(m_minimumVolumeSize));
  cout << "-------------------" << endl;
  cout << "minimumVolumeSize =" << m_minimumVolumeSize << endl;
  cout << "maxDistance =" << m_fNucleusSphereDiameter << endl;
  cout << "maxDepth = " <<  m_maxDepth << endl;
  cout << "mdr =" << m_dr << endl;
  
//  m_pOctree->addPoint(1146.27, 1325.35, 2017.46, 1, m_maxDepth);
//  x_ = 1145.66 ;
//  y_ = 1322.28 ;
//  z_ = 2092.46 ;


//  m_pOctree->addPoint(875.523,	-271.459,	608.49, 1, m_maxDepth);
//  x_ = 886.687;
//  y_ = -249.995;
//  z_ =	645.99;

  m_pOctree->addPoint(1256.19,  1478.53, -1194.18,1 , m_maxDepth);
  x_ = 1214.87;
  y_ = 1459.49;
  z_ = -1198.86;

//  m_pOctree->addPoint(-216.892    , 1984.53, 870.215,1,m_maxDepth);
//  x_ = -212.967;
//  y_ = 2036.59;
//  z_ = 870.215;

//  m_pOctree->addPoint(-578.7,  -2595.21,    -961.974,1, m_maxDepth );
//  x_ = -599.712;      
//  y_ = -2613.57; 
//  z_ = -957.286;

//  m_pOctree->addPoint(236.693, 995.11, 482.206, 1, m_maxDepth);
//  m_pOctree->addPoint(516.439, 319.745, -42.7936, 1, m_maxDepth);
//  m_pOctree->addPoint(196.918, 828.26, -713.106, 1, m_maxDepth);
//  m_pOctree->addPoint(356.531, 25.827, -1088.11, 1, m_maxDepth);
//  m_pOctree->addPoint(1104.85, 525.84, -1088.11, 1, m_maxDepth);
//  m_pOctree->addPoint(1829.41, 1009.98, -863.106, 1, m_maxDepth);
//  m_pOctree->addPoint(1972.03, 1726.94, -1388.11, 1, m_maxDepth);
//
//  m_pOctree->addPoint(1420.89, 2278.07, -938.106, 1, m_maxDepth);
//
//  m_pOctree->addPoint(1420.89, 1406.65, -713.106, 1, m_maxDepth);
//  m_pOctree->addPoint(996.227, 981.988, -1383.42, 1, m_maxDepth);
//  m_pOctree->addPoint(320.862, 702.243, -858.419, 1, m_maxDepth);
//  m_pOctree->addPoint(-12.6165, 1507.33, -1083.42, 1, m_maxDepth);
//  m_pOctree->addPoint(711.944, 1991.47, -1308.42, 1, m_maxDepth);
//  x_ = 1467.82;
//  y_ = 2304.56;
//  z_ = -933.419;

   minPosition_ = Eigen::Matrix<float, 4, 1> ( x_ - m_dr, y_ - m_dr, z_ - m_dr, 1);
   maxPosition_ = Eigen::Matrix<float, 4, 1> ( x_ + m_dr, y_ + m_dr, z_ + m_dr, 1);
  nodes_.clear();
  m_pOctree->getAllNodesInVolumeOnDepth(nodes_, minPosition_, maxPosition_, m_maxDepth, true);
  
  if (nodes_.size() > 0){
    bool bFlag = false;
    for (int i = 0; i< (int)nodes_.size(); i++) {
      Eigen::Matrix<float,4,1> point_ = (nodes_[i])->getPosition();
      float point_x_ = point_[0];
      float point_y_ = point_[1];
      float point_z_ = point_[2];
      float dist_ = sqrt((point_x_ - x_)*(point_x_ - x_) + (point_y_ - y_ )*(point_y_ - y_) + (point_z_ - z_)*(point_z_ - z_));
      if (dist_ < m_fCollisionLength) {
        bFlag = true;
        cout << x_ << " " << y_ << " " << z_ <<endl;
        cout << point_x_ << " " << point_y_ << " " << point_z_ <<endl;
        cout << dist_ << endl;
        cout << "collision happend";
        break;
      }
    }
    return bFlag;
  }else{
    cout << "aa" << endl;
    return false;
  }

  return 0;
  
  
}
