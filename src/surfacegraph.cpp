/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "surfacegraph.h"

#include <QtCore/qmath.h>

#include <QtDataVisualization/Q3DTheme>
#include <QtDataVisualization/QValue3DAxis>
#include <QtGui/QImage>

using namespace QtDataVisualization;

SurfaceGraph::SurfaceGraph(Q3DSurface *surface) : m_graph(surface) {
  m_graph->setAxisX(new QValue3DAxis);
  m_graph->setAxisY(new QValue3DAxis);
  m_graph->setAxisZ(new QValue3DAxis);

  //! [0]
  m_sqrtSinProxy = new QSurfaceDataProxy();
  m_sqrtSinSeries = new QSurface3DSeries(m_sqrtSinProxy);
  data_array_ = new QSurfaceDataArray();
  enableSqrtSinModel();
}

SurfaceGraph::~SurfaceGraph() { delete m_graph; }

void SurfaceGraph::Fill(const std::vector<double> &dd, size_t dim) {
  if (data_array_->size() != dim) {
    for (auto *item : *data_array_) {
      delete item;
    }
    data_array_->clear();
    data_array_->reserve(dim);
    for (size_t i = 0; i < dim; ++i) {
      *data_array_ << new QSurfaceDataRow(dim);
    }
  }

  size_t ind = 0;
  ymin_ = 99e9;
  ymax_ = -99e9;
  for (size_t i = 0; i < dim; ++i) {
    float x = static_cast<float>(i) / dim;
    QSurfaceDataRow *ptr_row_data = (*data_array_)[i];
    for (size_t j = 0; j < dim; ++j) {
      float z = static_cast<float>(j) / dim;
      float y = dd[ind++];
      (*ptr_row_data)[j].setPosition(QVector3D(x, y, z));

      ymin_ = std::min(y, ymin_);
      ymax_ = std::max(y, ymax_);
    }
  }
  m_sqrtSinProxy->resetArray(data_array_);
  m_graph->axisY()->setRange(ymin_, 10 * ymax_);
}

void SurfaceGraph::Adjust() {
  const float delta = 1.0 / data_array_->size();
  m_graph->axisX()->setRange(0.0f - delta, 1.0f + delta);
  m_graph->axisZ()->setRange(0.0f - delta, 1.0f + delta);
  m_graph->axisY()->setRange(ymin_, ymax_);
}

void SurfaceGraph::enableSqrtSinModel() {
  m_sqrtSinSeries->setDrawMode(QSurface3DSeries::DrawSurface);
  m_graph->addSeries(m_sqrtSinSeries);
  m_graph->axisX()->setRange(0.0f - 0.1, 1.1f);
  m_graph->axisZ()->setRange(0.0f - 0.1, 1.1f);
}

void SurfaceGraph::setAxisXRange(float min, float max) {
  m_graph->axisX()->setRange(min, max);
}

void SurfaceGraph::setAxisZRange(float min, float max) {
  m_graph->axisZ()->setRange(min, max);
}
//! [5]

//! [6]
void SurfaceGraph::changeTheme(int theme) {
  m_graph->activeTheme()->setType(Q3DTheme::Theme(theme));
}
//! [6]

void SurfaceGraph::setBlackToYellowGradient() {
  //! [7]
  QLinearGradient gr;
  gr.setColorAt(0.0, Qt::black);
  gr.setColorAt(0.33, Qt::blue);
  gr.setColorAt(0.67, Qt::red);
  gr.setColorAt(1.0, Qt::yellow);

  m_graph->seriesList().at(0)->setBaseGradient(gr);
  m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
  //! [7]
}

void SurfaceGraph::setGreenToRedGradient() {
  QLinearGradient gr;
  gr.setColorAt(0.0, Qt::darkGreen);
  gr.setColorAt(0.5, Qt::yellow);
  gr.setColorAt(0.8, Qt::red);
  gr.setColorAt(1.0, Qt::darkRed);

  m_graph->seriesList().at(0)->setBaseGradient(gr);
  m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
}

void SurfaceGraph::SetGradient(const QLinearGradient &gd) {
  m_graph->seriesList().at(0)->setBaseGradient(gd);
  // m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
  m_graph->seriesList().at(0)->setColorStyle(
      Q3DTheme::ColorStyleObjectGradient);
}
