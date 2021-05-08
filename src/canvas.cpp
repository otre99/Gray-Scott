#include "canvas.h"

#include <QPainter>

Canvas::Canvas(QWidget *parent) : QWidget(parent) { m_radius = 20; }

void Canvas::paintEvent(QPaintEvent * /*event*/) {
  QPainter p(this);
  p.drawImage(p.viewport(), m_dp_img);
}

void Canvas::mouseMoveEvent(QMouseEvent *event) { draw(event->pos()); }

void Canvas::mousePressEvent(QMouseEvent *event) {
  m_last_pt = event->pos();
  draw(m_last_pt);
}

void Canvas::mouseReleaseEvent(QMouseEvent *event) {}

void Canvas::clear() {
  m_dp_img = QImage(this->size(), QImage::Format_ARGB32);
  m_dp_img.fill(Qt::white);
  update();
}

void Canvas::enableEraser(bool enable) {
  if (enable) {
    m_draw_color = QColor(Qt::white);
  } else {
    m_draw_color = QColor(Qt::black);
  }
}

void Canvas::drawRandomCircles(int n) {
  if (m_dp_img.isNull()) return;
  QPainter p(&m_dp_img);
  p.setBrush(QBrush(Qt::black));
  for (int i = 0; i < n; ++i) {
    int x = std::rand() % m_dp_img.width();
    int y = std::rand() % m_dp_img.width();
    int rx = 1 + std::rand() % m_radius;
    int ry = 1 + std::rand() % m_radius;
    p.drawEllipse({x, y}, rx, ry);
  }
  update();
}

void Canvas::drawRandomSquares(int n) {
  if (m_dp_img.isNull()) return;
  QPainter p(&m_dp_img);
  p.setBrush(QBrush(Qt::black));
  for (int i = 0; i < n; ++i) {
    int x = std::rand() % m_dp_img.width();
    int y = std::rand() % m_dp_img.width();
    int rx = 1 + std::rand() % m_radius;
    int ry = 1 + std::rand() % m_radius;
    p.drawRect(x, y, rx, ry);
  }
  update();
}

void Canvas::draw(const QPoint &pt) {
  if (m_dp_img.size() != size()) {
    m_dp_img = m_dp_img.scaled(size());
  }

  QPainter p(&m_dp_img);
  p.setRenderHint(QPainter::Antialiasing);

  QPen pen;
  pen.setWidth(m_radius);
  pen.setColor(m_draw_color);
  p.setPen(pen);

  if (pt != m_last_pt) {
    p.drawLine(m_last_pt, pt);
  } else {
    p.drawPoint(pt);
  }
  m_last_pt = pt;
  update();
}
