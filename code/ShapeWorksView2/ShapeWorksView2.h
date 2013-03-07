#include <QMainWindow>

#include "tinyxml.h"

#include "itkParticleShapeLinearRegressionMatrixAttribute.h"
#include "itkParticlePositionReader.h"
#include "itkParticleShapeStatistics.h"
#include "itkParticlePositionWriter.h"

#include <vtkSmartPointer.h>

#include <ModelCache.h>

class vtkRenderer;
class vtkLookupTable;
class vtkColorTransferFunction;
class vtkPoints;
class vtkPolyData;
class vtkGlyph3D;
class vtkPolyDataMapper;
class vtkActor;
class vtkSphereSource;
class vtkUnsignedLongArray;
class vtkArrowSource;
class vtkTransform;

class CustomSurfaceReconstructionFilter;
class vtkContourFilter;
class vtkReverseSense;
class vtkSmoothPolyDataFilter;
class vtkPolyDataNormals;
class vtkDecimatePro;
class vtkImageConstantPad;
class vtkPowerCrustSurfaceReconstruction;

// Forward Qt class declarations
class Ui_ShapeWorksView2;

class ViewerLayout;

class ShapeWorksView2 : public QMainWindow
{
  Q_OBJECT
public:

  ShapeWorksView2( int argc, char** argv );
  ~ShapeWorksView2();

  virtual void closeEvent( QCloseEvent* event );

public Q_SLOTS:

  void on_actionQuit_triggered();
  void on_actionPreferences_triggered();

  // display modes
  void on_meanButton_clicked();
  void on_sampleButton_clicked();
  void on_pcaButton_clicked();

  // mean mode
  void on_meanOverallButton_clicked();
  void on_meanGroup1Button_clicked();
  void on_meanGroup2Button_clicked();

  // sample mode
  void on_sampleSpinBox_valueChanged();
  void on_medianButton_clicked();
  void on_medianGroup1Button_clicked();
  void on_medianGroup2Button_clicked();

  // PCA mode
  void on_pcaSlider_valueChanged();
  void on_pcaModeSpinBox_valueChanged();
  void on_pcaGroupSlider_valueChanged();

  // visualization
  void on_showGlyphs_stateChanged();
  void on_showSurface_stateChanged();
  void on_usePowerCrustCheckBox_stateChanged();
  void on_neighborhoodSpinBox_valueChanged();
  void on_spacingSpinBox_valueChanged();

private:

  void initializeRenderer();
  void initializeGlyphs();
  void initializeSurface();

  void updateShapeMode();
  void updateSurfaceSettings();
  void updateActors();

  void redraw();

  bool readParameterFile( char* filename );
  void displayShape( const vnl_vector<double> &pos );
  void computeModeShape();

  // designer form
  Ui_ShapeWorksView2* ui;

  vtkSmartPointer<vtkRenderer>             renderer;
  vtkSmartPointer<vtkLookupTable>          lut;
  vtkSmartPointer<vtkPoints>               glyphPoints;
  vtkSmartPointer<vtkPolyData>             glyphPointSet;
  vtkSmartPointer<vtkGlyph3D>              glyphs;
  vtkSmartPointer<vtkPolyDataMapper>       glyphMapper;
  vtkSmartPointer<vtkActor>                glyphActor;
  vtkSmartPointer<vtkSphereSource>         sphereSource;
  vtkSmartPointer<vtkUnsignedLongArray>    scalars;
  vtkSmartPointer<vtkPolyDataMapper>       surfaceMapper;
  vtkSmartPointer<vtkActor>                surfaceActor;
  vtkSmartPointer<vtkContourFilter>        surfaceContourFilter;
  vtkSmartPointer<vtkReverseSense>         surfaceReverseSense;
  vtkSmartPointer<vtkSmoothPolyDataFilter> surfaceSmoothFilter;
  vtkSmartPointer<vtkPolyDataNormals>      polydataNormals;

  vtkSmartPointer<CustomSurfaceReconstructionFilter>  surface;
  vtkSmartPointer<vtkPowerCrustSurfaceReconstruction> powercrust;

  //vtkSmartPointer<vtkColorTransferFunction>   differenceLUT;
  //vtkSmartPointer<vtkColorTransferFunction>   pValueTFunc;
  //vtkSmartPointer<vtkTransformPolyDataFilter> arrowFlipFilter;
  //vtkSmartPointer<vtkGlyph3D>                 arrowGlyphs;
  //vtkSmartPointer<vtkPolyDataNormals>         m_surfNormals;
  //vtkSmartPointer<vtkDecimatePro>             m_surfDecimate;
  //vtkSmartPointer<vtkTransform>               transform180;
  //vtkSmartPointer<vtkArrowSource>             arrowSource;

  ParticleShapeStatistics<3> stats;

  int numSamples;
  bool groupsAvailable;

  // a copy of the current shape mesh
  vnl_vector<double> currentShape;

  // cache of shape models
  ModelCache modelCache;
};