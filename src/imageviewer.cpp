#include "imageviewer.h"
#include <QDebug>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>
#include <QScrollBar>

ImageViewer::ImageViewer(QWidget *parent) : QAbstractScrollArea(parent) {
  Init();
  connect(horizontalScrollBar(), SIGNAL(valueChanged(int)), SLOT(SetXmov(int)));
  connect(verticalScrollBar(), SIGNAL(valueChanged(int)), SLOT(SetYmov(int)));
  viewport()->setMouseTracking(true);
  viewport()->setCursor(QCursor(Qt::CrossCursor));
}

void ImageViewer::AttachImagePtr(QImage *ptr) {
  image_ptr_ = ptr;
  const double s1 = double(viewport()->width()) / image_ptr_->width();
  const double s2 = double(viewport()->height()) / image_ptr_->height();
  scf_ = 1.0; //std::min(std::min(s1, s2), 1.0);
  AdjustAll();
}

void ImageViewer::Init() {
  image_ptr_ = nullptr;
  scf_ = 1.0;
  xmov_ = ymov_ = 0;
}

void ImageViewer::resizeEvent(QResizeEvent *) {
  if (image_ptr_)
    AdjustAll();
}

void ImageViewer::paintEvent(QPaintEvent *) {
  if (!image_ptr_)
    return;
  QPainter p(viewport());
  p.translate(QPoint(viewport()->width() / 2, viewport()->height() / 2));
  p.drawImage(QRect(-screen_w_ / 2, -screen_h_ / 2, screen_w_, screen_h_),
              *image_ptr_, QRect(xmov_, ymov_, cw_, ch_));

}

void ImageViewer::mouseMoveEvent(QMouseEvent *e) {
  const QPoint pd_pos = viewport()->mapFrom(this, e->pos());
  const int xx =
      xmov_ + (pd_pos.x() - viewport()->width() / 2 + screen_w_ / 2) / scf_;
  const int yy =
      ymov_ + (pd_pos.y() - viewport()->height() / 2 + screen_h_ / 2) / scf_;
  emit PixelTrack(xx, yy);
}

void ImageViewer::AdjustAll() {
  screen_w_ = viewport()->width();
  screen_h_ = viewport()->height();
  cw_ = screen_w_ / scf_ + .5;
  if (cw_ > image_ptr_->width()) {
    cw_ = image_ptr_->width();
    screen_w_ = cw_ * scf_ + .5;
  }
  ch_ = screen_h_ / scf_ + .5;
  if (ch_ > image_ptr_->height()) {
    ch_ = image_ptr_->height();
    screen_h_ = ch_ * scf_ + .5;
  }
  horizontalScrollBar()->setPageStep(cw_);
  horizontalScrollBar()->setMaximum(image_ptr_->width() - cw_);

  verticalScrollBar()->setPageStep(ch_);
  verticalScrollBar()->setMaximum(image_ptr_->height() - ch_);
  viewport()->update();
}

void ImageViewer::FixWidth() {
  SetScf(double(viewport()->width()) / image_ptr_->width());
}

void ImageViewer::SetXmov(int x) {
  if (x != xmov_) {
    xmov_ = x;
    viewport()->update();
  }
}

void ImageViewer::SetYmov(int y) {
  if (y != ymov_) {
    ymov_ = y;
    viewport()->update();
  }
}

void ImageViewer::SetScf(const double s) {
  if (s != scf_) {
    scf_ = s;
    AdjustAll();
  }
}

QImage *ImageViewer::ImagePtr() { return image_ptr_; }
