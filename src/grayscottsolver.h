#ifndef GRAYSCOTTSOLVER_H
#define GRAYSCOTTSOLVER_H
#include <QThread>
#include <vector>

class GrayScottSolver : public QObject {
  Q_OBJECT
  QThread *thread_;

 public:
  GrayScottSolver();
  ~GrayScottSolver() override;
  void initParams(int dim, double du, double dv, double f, double k);
  void initialize();
  void initializeRandom();
  void initialize(const QImage &img);
  void updateDpData();
  size_t getSize() { return m_dim; }
  const std::vector<double> &getData() const { return m_data_dp; }
  bool isRunning() { return !loop_break_; }
  void selectUOrV(bool select_u) { m_u_is_selected = select_u; }
  void selectMethod(int method) { m_method = method; }
  void EnableParallel(bool enable) { enable_parallel = enable; }

 public slots:
  void stop() { loop_break_ = true; }
  void solve(int n, double dt);

 private:
  void fEuler(int n);
  void eval(const double *u, const double *v, double *ou, double *ov,
            size_t size);
  void cNSplitting(int n);
  double sourceU(double u, double v);
  double sourceV(double u, double v);

  std::vector<double> m_u;
  std::vector<double> m_v;
  std::vector<double> m_k1u, m_k1v;
  std::vector<double> m_data_dp;
  bool loop_break_;
  double m_f, m_k;
  double m_du, m_dv;
  double m_dt;
  size_t m_dim;
  bool m_u_is_selected;
  int m_method;
  bool enable_parallel;

 signals:
  void dataReady();
};

#endif  // GRAYSCOTTSOLVER_H
