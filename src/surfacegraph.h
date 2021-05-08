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

#ifndef SURFACEGRAPH_H
#define SURFACEGRAPH_H

#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QHeightMapSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtWidgets/QSlider>

using namespace QtDataVisualization;

class SurfaceGraph : public QObject {
  Q_OBJECT
 public:
  explicit SurfaceGraph(Q3DSurface *surface);
  ~SurfaceGraph();

  void enableSqrtSinModel();

  //! [0]
  void toggleModeNone() {
    m_graph->setSelectionMode(QAbstract3DGraph::SelectionNone);
  }
  void toggleModeItem() {
    m_graph->setSelectionMode(QAbstract3DGraph::SelectionItem);
  }
  void toggleModeSliceRow() {
    m_graph->setSelectionMode(QAbstract3DGraph::SelectionItemAndRow |
                              QAbstract3DGraph::SelectionSlice);
  }
  void toggleModeSliceColumn() {
    m_graph->setSelectionMode(QAbstract3DGraph::SelectionItemAndColumn |
                              QAbstract3DGraph::SelectionSlice);
  }
  //! [0]

  void setBlackToYellowGradient();
  void setGreenToRedGradient();
  void SetGradient(const QLinearGradient &gd);
  void Fill(const std::vector<double> &dd, size_t dim);
  void Adjust();

 public Q_SLOTS:
  void changeTheme(int theme);

 private:
  QSurfaceDataArray *data_array_;
  Q3DSurface *m_graph;
  QSurfaceDataProxy *m_sqrtSinProxy;
  QSurface3DSeries *m_sqrtSinSeries;
  void setAxisXRange(float min, float max);
  void setAxisZRange(float min, float max);
  float ymin_, ymax_;
};

#endif  // SURFACEGRAPH_H
