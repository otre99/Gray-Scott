#ifndef CANVAS_H
#define CANVAS_H
#include <QMouseEvent>
#include <QPaintEvent>
#include <QWidget>

class Canvas : public QWidget {
  Q_OBJECT
 public:
  Canvas(QWidget *parent = nullptr);
  void paintEvent(QPaintEvent *event) override;
  void mouseMoveEvent(QMouseEvent *event) override;
  void mousePressEvent(QMouseEvent *event) override;
  void mouseReleaseEvent(QMouseEvent *event) override;
  void setRadius(int radius) { m_radius = radius; }
  void enableEraser(bool enable);
  QImage &getImage() { return m_dp_img; }

 public slots:
  void clear();
  void drawRandomCircles(int n);
  void drawRandomSquares(int n);

 private:
  void draw(const QPoint &pt);
  QImage m_dp_img;
  QPoint m_last_pt;
  QColor m_draw_color;
  int m_radius;
};

#endif  // CANVAS_H
