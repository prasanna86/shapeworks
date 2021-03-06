// generated by Fast Light User Interface Designer (fluid) version 1.0110

#ifndef AnalyzeCorrespondenceGUI_h
#define AnalyzeCorrespondenceGUI_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include "vtkFlRenderWindowInteractor.h"
#include <FL/Fl_Group.H>
#include <FL/Fl_Spinner.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Value_Output.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Box.H>

class AnalyzeCorrespondenceGUI {
public:
  AnalyzeCorrespondenceGUI();
  Fl_Double_Window *mainWindow;
  Fl_Menu_Bar *menu;
  static Fl_Menu_Item menu_menu[];
private:
  void cb_Load_i(Fl_Menu_*, void*);
  static void cb_Load(Fl_Menu_*, void*);
  void cb_Load1_i(Fl_Menu_*, void*);
  static void cb_Load1(Fl_Menu_*, void*);
  void cb_Load2_i(Fl_Menu_*, void*);
  static void cb_Load2(Fl_Menu_*, void*);
  void cb_Load3_i(Fl_Menu_*, void*);
  static void cb_Load3(Fl_Menu_*, void*);
  void cb_Write_i(Fl_Menu_*, void*);
  static void cb_Write(Fl_Menu_*, void*);
  void cb_Write1_i(Fl_Menu_*, void*);
  static void cb_Write1(Fl_Menu_*, void*);
  void cb_Point_i(Fl_Menu_*, void*);
  static void cb_Point(Fl_Menu_*, void*);
  void cb_Exit_i(Fl_Menu_*, void*);
  static void cb_Exit(Fl_Menu_*, void*);
  void cb_Samples_i(Fl_Menu_*, void*);
  static void cb_Samples(Fl_Menu_*, void*);
  void cb_Mean_i(Fl_Menu_*, void*);
  static void cb_Mean(Fl_Menu_*, void*);
  void cb_Median_i(Fl_Menu_*, void*);
  static void cb_Median(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *menuGroup1Mean;
private:
  void cb_menuGroup1Mean_i(Fl_Menu_*, void*);
  static void cb_menuGroup1Mean(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *menuGroup1Median;
private:
  void cb_menuGroup1Median_i(Fl_Menu_*, void*);
  static void cb_menuGroup1Median(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *menuGroup2Mean;
private:
  void cb_menuGroup2Mean_i(Fl_Menu_*, void*);
  static void cb_menuGroup2Mean(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *menuGroup2Median;
private:
  void cb_menuGroup2Median_i(Fl_Menu_*, void*);
  static void cb_menuGroup2Median(Fl_Menu_*, void*);
  void cb_PCA_i(Fl_Menu_*, void*);
  static void cb_PCA(Fl_Menu_*, void*);
public:
  static Fl_Menu_Item *ViewLinearRegressionChoice;
private:
  void cb_ViewLinearRegressionChoice_i(Fl_Menu_*, void*);
  static void cb_ViewLinearRegressionChoice(Fl_Menu_*, void*);
  void cb_Mean1_i(Fl_Menu_*, void*);
  static void cb_Mean1(Fl_Menu_*, void*);
  void cb_Reconstruction_i(Fl_Menu_*, void*);
  static void cb_Reconstruction(Fl_Menu_*, void*);
  void cb_Color_i(Fl_Menu_*, void*);
  static void cb_Color(Fl_Menu_*, void*);
  void cb_Glyphs_i(Fl_Menu_*, void*);
  static void cb_Glyphs(Fl_Menu_*, void*);
public:
  vtkFlRenderWindowInteractor *imageView;
  Fl_Group *OptionsReconstructionGroup;
  Fl_Spinner *samplespacing;
private:
  void cb_samplespacing_i(Fl_Spinner*, void*);
  static void cb_samplespacing(Fl_Spinner*, void*);
public:
  Fl_Spinner *neighborhoodsize;
private:
  void cb_neighborhoodsize_i(Fl_Spinner*, void*);
  static void cb_neighborhoodsize(Fl_Spinner*, void*);
  void cb_ON_i(Fl_Button*, void*);
  static void cb_ON(Fl_Button*, void*);
  void cb_OFF_i(Fl_Button*, void*);
  static void cb_OFF(Fl_Button*, void*);
public:
  Fl_Spinner *smooth_iter;
private:
  void cb_smooth_iter_i(Fl_Spinner*, void*);
  static void cb_smooth_iter(Fl_Spinner*, void*);
public:
  Fl_Spinner *gradsmooth;
  Fl_Group *OptionsColorSchemeGroup;
  Fl_Spinner *color_switcher;
private:
  void cb_color_switcher_i(Fl_Spinner*, void*);
  static void cb_color_switcher(Fl_Spinner*, void*);
public:
  Fl_Group *OptionsGlyphGroup;
  Fl_Value_Slider *glyph_scale;
private:
  void cb_glyph_scale_i(Fl_Value_Slider*, void*);
  static void cb_glyph_scale(Fl_Value_Slider*, void*);
public:
  Fl_Value_Slider *glyph_quality;
private:
  void cb_glyph_quality_i(Fl_Value_Slider*, void*);
  static void cb_glyph_quality(Fl_Value_Slider*, void*);
public:
  Fl_Light_Button *show_correspondence_button;
private:
  void cb_show_correspondence_button_i(Fl_Light_Button*, void*);
  static void cb_show_correspondence_button(Fl_Light_Button*, void*);
  void cb_ON1_i(Fl_Button*, void*);
  static void cb_ON1(Fl_Button*, void*);
  void cb_OFF1_i(Fl_Button*, void*);
  static void cb_OFF1(Fl_Button*, void*);
public:
  Fl_Group *ViewMeanGroup;
  Fl_Group *ViewPCAModesGroup;
  Fl_Value_Slider *mode_position;
private:
  void cb_mode_position_i(Fl_Value_Slider*, void*);
  static void cb_mode_position(Fl_Value_Slider*, void*);
public:
  Fl_Value_Output *lambda_display;
  Fl_Slider *groupdiff_position;
private:
  void cb_groupdiff_position_i(Fl_Slider*, void*);
  static void cb_groupdiff_position(Fl_Slider*, void*);
public:
  Fl_Box *group1LabelBox;
  Fl_Box *group2LabelBox;
  Fl_Spinner *mode;
private:
  void cb_mode_i(Fl_Spinner*, void*);
  static void cb_mode(Fl_Spinner*, void*);
public:
  Fl_Value_Slider *simple_regression;
private:
  void cb_simple_regression_i(Fl_Value_Slider*, void*);
  static void cb_simple_regression(Fl_Value_Slider*, void*);
public:
  Fl_Value_Output *loading_display;
  Fl_Group *ViewRegressionGroup;
  Fl_Value_Slider *position;
private:
  void cb_position_i(Fl_Value_Slider*, void*);
  static void cb_position(Fl_Value_Slider*, void*);
public:
  Fl_Group *ViewSampleGroup;
  Fl_Spinner *sample_selector;
private:
  void cb_sample_selector_i(Fl_Spinner*, void*);
  static void cb_sample_selector(Fl_Spinner*, void*);
public:
  Fl_Value_Output *group_display;
  virtual ~AnalyzeCorrespondenceGUI();
  virtual void Quit();
  virtual void Show();
  virtual void Hide();
  virtual void ComputeRegressionShape();
  virtual void SetGlyphScale();
  virtual void ComputeSurface();
  virtual void RemoveSurface();
  virtual void ChangeColorScheme();
  virtual void HideGroups();
  virtual void DisplayMean(unsigned int);
  virtual void DisplaySamples();
  virtual void ComputeModeShape();
  virtual void ComputeGroupMeanDifferenceShape();
  virtual void WritePoints();
  virtual void WritePCALoadings();
  virtual void LoadScalars();
  virtual void DisplayMeanDifference();
  virtual void DisplayGroupMedian(int);
  virtual void ShowCorrespondence();
  virtual void ComputeSimpleRegressionShape();
  virtual void ComputeSimpleRegressionParameters();
  virtual void LoadPCAShape();
  virtual void LoadPointFile();
  virtual void ShowSpheres();
  virtual void SetSurfaceSmoothing();
  virtual void LoadVectorField();
  virtual void PointFileDiff();
  virtual void ShowGlyphs();
  virtual void RemoveGlyphs();
};
#endif
