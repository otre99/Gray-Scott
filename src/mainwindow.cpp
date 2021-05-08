#include "mainwindow.h"

#include <QDebug>
#include <QMap>
#include <QPainter>
#include <cmath>

#include "./ui_mainwindow.h"

static QMap<QString, QCPColorGradient::GradientPreset> name2gp;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  SetUp();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_push_button_run__clicked() {
  if (!m_gray_scott_solver.isRunning()) {
    emit StartSolver(ui->spin_box_iters_->value(),
                     ui->line_dt_->text().toDouble());
    ui->push_button_run_->setText("Stop");
    ui->group_box_params_->setEnabled(false);
    ui->groupBox_3->setEnabled(false);
  } else {
    m_gray_scott_solver.stop();
    ui->push_button_run_->setText("Run");
    ui->group_box_params_->setEnabled(true);
    ui->groupBox_3->setEnabled(true);
  }
}

void MainWindow::on_push_button_init__clicked() {
  UpdateParams();
  m_gray_scott_solver.initialize();
  UpdatePlot();
}

void MainWindow::on_groupBox_3_clicked(bool checked) {
  if (checked) {
    setCentralWidget(m_canvas = new Canvas);
    m_canvas->clear();
    m_canvas->setRadius(ui->spin_box_radius_->value());
    ui->group_box_plot_->setEnabled(false);
    ui->group_box_params_->setEnabled(false);
    ui->group_box_solver_->setEnabled(false);
  } else {
    UpdateParams();
    m_gray_scott_solver.initialize(m_canvas->getImage());
    if (m_active_dp == kDp2D) {
      m_active_dp = kDp3D;
    } else {
      m_active_dp = kDp2D;
    }
    SwitchDP();
    // UpdatePlot();
    ui->group_box_plot_->setEnabled(true);
    ui->group_box_params_->setEnabled(true);
    ui->group_box_solver_->setEnabled(true);
  }
}

void MainWindow::on_spin_box_radius__valueChanged(int arg1) {
  m_canvas->setRadius(arg1);
}

void MainWindow::on_pb_adjust__clicked() {
  if (m_active_dp == kDp2D) {
    m_custom_plot->rescaleAxes();
    m_custom_plot->replot();
  } else {
    m_surface_modifier->Adjust();
  }
}

void MainWindow::on_pb_rescale__clicked() {
  if (m_active_dp == kDp2D) {
    m_color_map->rescaleDataRange();
    m_custom_plot->replot();
  }
}

void MainWindow::on_combo_box_cmaps__activated(const QString &arg1) {
  if (m_active_dp == kDp3D) {
    const QCPColorGradient qcgd(name2gp[arg1]);
    QLinearGradient gd;
    const auto cv = qcgd.colorStops();
    for (const auto key : cv.keys()) {
      gd.setColorAt(key, cv[key]);
    }
    m_surface_modifier->SetGradient(gd);
  } else {
    m_color_map->setGradient(name2gp[arg1]);
    m_custom_plot->replot();
  }
}

void MainWindow::on_check_box_log__clicked(bool checked) {
  if (m_active_dp == kDp2D) {
    if (checked) {
      m_color_map->colorScale()->setDataScaleType(QCPAxis::stLogarithmic);
    } else {
      m_color_map->colorScale()->setDataScaleType(QCPAxis::stLinear);
    }
    m_custom_plot->replot();
  }
}

void MainWindow::on_tool_button_swith_dp__clicked() {
  SwitchDP();
  UpdatePlot();
  m_color_map->rescaleDataRange();
}

void MainWindow::on_combo_box_uov__activated(int index) {
  m_gray_scott_solver.selectUOrV(!index);
  if (!m_gray_scott_solver.isRunning()) {
    m_gray_scott_solver.updateDpData();
    UpdatePlot();
  }
}

void MainWindow::on_combo_box_method__activated(int index) {
  m_gray_scott_solver.selectMethod(index);
}

void MainWindow::on_check_box_eraser__clicked(bool checked) {
  m_canvas->enableEraser(checked);
}

void MainWindow::on_push_button_circles__clicked() {
  m_canvas->setRadius(ui->spin_box_radius_->value());
  m_canvas->drawRandomCircles(ui->spin_box_n_->value());
}

void MainWindow::on_push_button_clear__clicked() { m_canvas->clear(); }

void MainWindow::on_push_button_squares__clicked() {
  m_canvas->setRadius(ui->spin_box_radius_->value());
  m_canvas->drawRandomSquares(ui->spin_box_n_->value());
}

void MainWindow::on_check_box_allow_parallel__clicked(bool checked) {
  m_gray_scott_solver.EnableParallel(checked);
}

void MainWindow::FillColormap(const std::vector<double> &dd, size_t dim) {
  m_color_map->data()->setSize(dim, dim);
  m_color_map->data()->setRange(QCPRange(0, 1), QCPRange(0, 1));
  size_t k = 0;
  for (int yIndex = dim - 1; yIndex >= 0; --yIndex) {
    for (size_t xIndex = 0; xIndex < dim; ++xIndex) {
      m_color_map->data()->setCell(xIndex, yIndex, dd[k++]);
    }
  }
  m_custom_plot->replot();
}

void MainWindow::SwitchDP() {
  if (m_active_dp == kDp3D) {
    m_active_dp = kDp2D;
    setCentralWidget(m_custom_plot = new QCustomPlot);
    SetupColorMapPlot(m_custom_plot);
  } else {
    m_surface_graph = new Q3DSurface();
    QWidget *container = QWidget::createWindowContainer(m_surface_graph);
    m_active_dp = kDp3D;
    m_surface_modifier = new SurfaceGraph(m_surface_graph);
    setCentralWidget(container);
  }
  on_combo_box_cmaps__activated(ui->combo_box_cmaps_->currentText());
  on_combo_box_uov__activated(ui->combo_box_uov_->currentIndex());
}

void MainWindow::SetupColorMapPlot(QCustomPlot *customPlot) {
  customPlot->setInteractions(
      QCP::iRangeDrag | QCP::iRangeZoom);  // this will also allow rescaling the
                                           // color scale by dragging/zooming
  customPlot->axisRect()->setupFullAxesBox(true);
  customPlot->xAxis->setLabel("x");
  customPlot->yAxis->setLabel("y");

  m_color_map = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
  m_color_map->data()->setRange(QCPRange(0, 1), QCPRange(0, 1));

  auto *colorScale = new QCPColorScale(customPlot);
  customPlot->plotLayout()->addElement(
      0, 1, colorScale);  // add it to the right of the main axis rect
  colorScale->setType(
      QCPAxis::atRight);  // scale shall be vertical bar with tick/axis labels
                          // right (actually atRight is already the default)
  m_color_map->setColorScale(
      colorScale);  // associate the color map with the color scale

  m_color_map->setGradient(QCPColorGradient::gpGrayscale);
  auto *marginGroup = new QCPMarginGroup(customPlot);
  customPlot->axisRect()->setMarginGroup(QCP::msBottom | QCP::msTop,
                                         marginGroup);
  colorScale->setMarginGroup(QCP::msBottom | QCP::msTop, marginGroup);

  on_pb_adjust__clicked();
}

void MainWindow::SetUp() {
  name2gp["gpGrayscale"] = QCPColorGradient::gpGrayscale;
  name2gp["gpHot"] = QCPColorGradient::gpHot;
  name2gp["gpCold"] = QCPColorGradient::gpCold;
  name2gp["gpNight"] = QCPColorGradient::gpNight;
  name2gp["gpCandy"] = QCPColorGradient::gpCandy;
  name2gp["gpGeography"] = QCPColorGradient::gpGeography;
  name2gp["gpIon"] = QCPColorGradient::gpIon;
  name2gp["gpThermal"] = QCPColorGradient::gpThermal;
  name2gp["gpPolar"] = QCPColorGradient::gpPolar;
  name2gp["gpSpectrum"] = QCPColorGradient::gpSpectrum;
  name2gp["gpJet"] = QCPColorGradient::gpJet;
  name2gp["gpHues"] = QCPColorGradient::gpHues;

  m_active_dp = kDp3D;
  SwitchDP();

  connect(&m_gray_scott_solver, &GrayScottSolver::dataReady, this,
          &MainWindow::UpdatePlot);
  connect(this, &MainWindow::StartSolver, &m_gray_scott_solver,
          &GrayScottSolver::solve);
  UpdateParams();
}

void MainWindow::UpdateParams() {
  double du = ui->le_du_->text().toDouble();
  double dv = ui->le_dv_->text().toDouble();
  double f = ui->le_f_->text().toDouble();
  double k = ui->le_k_->text().toDouble();
  m_gray_scott_solver.initParams(ui->le_res_->text().toUInt(), du, dv, f, k);
}

void MainWindow::UpdatePlot() {
  if (m_active_dp == kDp2D) {
    FillColormap(m_gray_scott_solver.getData(), m_gray_scott_solver.getSize());
  } else {
    m_surface_modifier->Fill(m_gray_scott_solver.getData(),
                             m_gray_scott_solver.getSize());
  }
}
