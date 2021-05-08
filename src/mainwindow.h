#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "canvas.h"
#include "grayscottsolver.h"
#include "qcustomplot.h"
#include "surfacegraph.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
  Q_OBJECT

  enum { kDp2D, kDp3D } m_active_dp = kDp2D;

 public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

 private slots:
  void on_push_button_run__clicked();
  void on_push_button_init__clicked();
  void on_groupBox_3_clicked(bool checked);
  void on_spin_box_radius__valueChanged(int arg1);
  void on_pb_adjust__clicked();
  void on_pb_rescale__clicked();
  void on_combo_box_cmaps__activated(const QString &arg1);
  void on_check_box_log__clicked(bool checked);
  void on_tool_button_swith_dp__clicked();
  void on_combo_box_uov__activated(int index);
  void on_combo_box_method__activated(int index);
  void on_check_box_eraser__clicked(bool checked);
  void on_push_button_circles__clicked();
  void on_push_button_clear__clicked();
  void on_push_button_squares__clicked();
  void on_check_box_allow_parallel__clicked(bool checked);

 private:
  // void GenColorMap(size_t res=1000);
  void FillColormap(const std::vector<double> &dd, size_t dim);
  Ui::MainWindow *ui;
  GrayScottSolver m_gray_scott_solver;
  Canvas *m_canvas;

  Q3DSurface *m_surface_graph;
  SurfaceGraph *m_surface_modifier;

  void SwitchDP();

  QCustomPlot *m_custom_plot;
  QCPColorMap *m_color_map;
  void SetupColorMapPlot(QCustomPlot *customPlot);
  void SetUp();
  void UpdateParams();

 signals:
  void StartSolver(int, double);
  void StopSolver();

 public slots:
  void UpdatePlot();
};
#endif  // MAINWINDOW_H
