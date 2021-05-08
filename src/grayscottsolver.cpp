#include "grayscottsolver.h"
#include <QApplication>
#include <QDebug>
#include <QImage>
#include <algorithm>
#include <future>
#include <random>

#include "linalg.h"



constexpr double delta = 0.05;


GrayScottSolver::GrayScottSolver()
    : loop_break_(true), m_dt(1.0), m_dim(0), m_u_is_selected(true), m_method(0) {
    thread_ = new QThread();
    this->moveToThread(thread_);
    thread_->start();
}

GrayScottSolver::~GrayScottSolver() {
    loop_break_ = true;
    thread_->quit();
    thread_->wait();
    delete thread_;
}

void GrayScottSolver::updateDpData() {
    if (m_u_is_selected)
        m_data_dp = m_u;
    else
        m_data_dp = m_v;
}


void GrayScottSolver::initParams(int dim, double du, double dv, double f,
                                 double k) {

    m_f = f;
    m_k = k;
    m_du = du;
    m_dv = dv;
    m_dim = dim;
    m_u.resize(dim * dim);
    m_v.resize(dim * dim);

    m_k1u.resize(dim * dim);
    m_k1v.resize(dim * dim);
}


void GrayScottSolver::initialize(const QImage &img) {
    initializeRandom();
    QImage dd = img.scaled(m_dim, m_dim);

    QRgb *img_data = reinterpret_cast<QRgb *>(dd.bits());
    QRgb white = QColor(Qt::white).rgb();

    for (size_t i = 0; i < m_dim * m_dim; ++i) {
        if (*img_data != white) {
            m_u[i] += 0.50;
            m_v[i] += 0.25;

        }
        ++img_data;
    }
    m_data_dp = m_u;
}

void GrayScottSolver::initializeRandom() {

    // TODO(rbt): allow more initializations
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, delta);

    std::generate(m_u.begin(), m_u.end(),
                  [&dist, &mt]() { return 1 + dist(mt); });
    std::generate(m_v.begin(), m_v.end(), [&dist, &mt]() { return dist(mt); });
}

void GrayScottSolver::initialize() {
    initializeRandom();

    const size_t radius = m_dim / 3;
    for (size_t i = radius; i < 2 * radius; ++i) {

        size_t i1 = m_dim * i + radius;
        size_t i2 = i1 + radius;

        std::fill(m_u.begin() + i1, m_u.begin() + i2, 0.5);
        std::fill(m_v.begin() + i1, m_v.begin() + i2, 0.25);
    }
    updateDpData();
}


void GrayScottSolver::solve(int n, double dt) {
    m_dt = dt;
    loop_break_ = false;
    while (true) {
        if (loop_break_) break;
        switch (m_method) {
        case 0:
            fEuler(n);
            break;
        case 1:
            cNSplitting(n);
            break;

        }
    }
}


void GrayScottSolver::fEuler(int n) {
    for (size_t i = 0; i < n; ++i) {
        eval(m_u.data(), m_v.data(), m_k1u.data(), m_k1v.data(), m_dim * m_dim);
        #pragma omp parallel for
        for (size_t j = 0; j < m_u.size(); ++j) {
            m_u[j] += m_dt * m_k1u[j];
            m_v[j] += m_dt * m_k1v[j];
        }
        if (loop_break_)
            break;
    }
    updateDpData();
    emit dataReady();
}


void SolveXDirection(const double *au, const double *bu, const double *cu,
                     const double *av, const double *bv, const double *cv, const double uf, double vf, const int N, double *u, double *v)
{
#pragma omp parallel for
    for (int k=0; k<N; ++k){
        double *pu = u+k*N; double *ru = new double [N];
        double *pv = v+k*N; double *rv = new double [N];
        int j=0;
        ru[j] = (pu[j+1]+pu[N-1])*uf + (1-2*uf)*pu[j];
        rv[j] = (pv[j+1]+pv[N-1])*vf + (1-2*vf)*pv[j];
        for (j=1; j<N-1; ++j) {
            ru[j] = (pu[j+1]+pu[j-1])*uf + (1-2*uf)*pu[j];
            rv[j] = (pv[j+1]+pv[j-1])*vf + (1-2*vf)*pv[j];
        }
        ru[j] = (pu[0]+pu[j-1])*uf + (1-2*uf)*pu[j];
        rv[j] = (pv[0]+pv[j-1])*vf + (1-2*vf)*pv[j];
        LinalAlg::CyclicPackV1(au,bu,cu,N,ru,pu);
        LinalAlg::CyclicPackV1(av,bv,cv,N,rv,pv);
        delete [] ru;
        delete [] rv;
    }
}

void SolveYDirection(const double *au, const double *bu, const double *cu,
                     const double *av, const double *bv, const double *cv, const double uf, double vf, const int N, double *u, double *v)
{
#pragma omp parallel for
    for (int k=0; k<N; ++k){
        double *pu = u+k;
        double *ru = new double [N];
        double *pv = v+k;
        double *rv = new double [N];
        const int NL=(N-1)*N;

        int j=0;
        int jj = j*N;
        ru[j] = (pu[jj+N]+pu[NL])*uf + (1-2*uf)*pu[jj];
        rv[j] = (pv[jj+N]+pv[NL])*vf + (1-2*vf)*pv[jj];
        for (j=1; j<N-1; ++j) {
            jj=j*N;
            ru[j] = (pu[jj+N]+pu[jj-N])*uf + (1-2*uf)*pu[jj];
            rv[j] = (pv[jj+N]+pv[jj-N])*vf + (1-2*vf)*pv[jj];
        }
        jj=j*N;
        ru[j] = (pu[jj-N]+pu[0])*uf + (1-2*uf)*pu[jj];
        rv[j] = (pv[jj-N]+pv[0])*vf + (1-2*vf)*pv[jj];

        LinalAlg::CyclicPackV1(au,bu,cu,N,ru,ru);
        LinalAlg::CyclicPackV1(av,bv,cv,N,rv,rv);
        for (j=0; j<N; ++j) {
            jj = j*N;
            pu[jj] = ru[j];
            pv[jj] = rv[j];
        }
        delete [] ru;
        delete [] rv;
    }
}



void GrayScottSolver::cNSplitting(int n) {
    const int N = std::sqrt(m_u.size());
    const double dt = 0.5*m_dt;
    double *au = new double [N];
    double *bu = new double [N];
    double *cu = new double [N];
    double *av = new double [N]; std::fill(au, au+N,m_du);
    double *bv = new double [N]; std::fill(bu, bu+N,m_du);
    double *cv = new double [N]; std::fill(cu, cu+N,m_du);
    const double uf = 0.5*dt*m_du;
    const double vf = 0.5*dt*m_dv;
    for (int i=0; i<N; ++i){
        au[i] = cu[i] = -uf; bu[i] = 1+2*uf;
        av[i] = cv[i] = -vf; bv[i] = 1+2*vf;
    }
    for (int i=0; i<n; ++i){
        // x-direction
        SolveXDirection(au,bu,cu,av,bv,cv,uf,vf, N, m_u.data(), m_v.data());
        // y-direction
        SolveYDirection(au,bu,cu,av,bv,cv,uf, vf,N,m_u.data(), m_v.data());
        // nonlinear part
        #pragma omp parallel for
        for (int k=0; k<N*N;++k){
            double uv2 = m_u[k]*m_v[k]*m_v[k];
            double uk1 = m_u[k] + dt*(-uv2 + m_f*(1-m_u[k]));
            double vk1 = m_v[k] + dt*(uv2 - (m_f+m_k)*m_v[k]);

            uv2 = uk1*vk1*vk1;
            double uk2 = -uv2 + m_f*(1-uk1);
            double vk2 =  uv2 - (m_f+m_k)*vk1;
            m_u[k]+=2*dt*uk2;
            m_v[k]+=2*dt*vk2;
        }
        // y-direction
        SolveYDirection(au,bu,cu,av,bv,cv,uf, vf,N,m_u.data(), m_v.data());
        // x-direction
        SolveXDirection(au,bu,cu,av,bv,cv,uf, vf,N,m_u.data(), m_v.data());

        if (loop_break_)
            break;
    }
    delete [] au; delete [] av;
    delete [] bu; delete [] bv;
    delete [] cu; delete [] cv;
    updateDpData();
    emit dataReady();
}

void GrayScottSolver::eval(const double *u, const double *v, double *ou,
                           double *ov, size_t size) {
    //#pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        double tu, tv, lu, lv, ru, rv, bu, bv;
        const int row = i / m_dim;
        const int col = i % m_dim;

        const int row_T = (row-1)%m_dim;
        const int row_B = (row+1)%m_dim;
        const int col_L = (col-1)%m_dim;
        const int col_R = (col+1)%m_dim;

        tu = u[row_T*m_dim+col];
        bu = u[row_B*m_dim+col];
        lu = u[col_L*m_dim+col];
        ru = u[col_R*m_dim+col];

        tv = v[row_T*m_dim+col];
        bv = v[row_B*m_dim+col];
        lv = v[col_L*m_dim+col];
        rv = v[col_R*m_dim+col];

        if (row) {
            tu = u[i - m_dim];
            tv = v[i - m_dim];
        } else {
            tu = tv = 0.0;
        }

        if (row != m_dim - 1) {
            bu = u[i + m_dim];
            bv = v[i + m_dim];
        } else {
            bu = bv = 0.0;
        }

        if (col) {
            lu = u[i - 1];
            lv = v[i - 1];
        } else {
            lu = lv = 0.0;
        }

        if (col != m_dim - 1) {
            ru = u[i + 1];
            rv = v[i + 1];
        } else {
            ru = rv = 0.0;
        }
        ou[i] = m_du * (-4 * u[i] + tu + bu + lu + ru) + sourceU(u[i], v[i]);
        ov[i] = m_dv * (-4 * v[i] + tv + bv + lv + rv) + sourceV(u[i], v[i]);
    }
}

double GrayScottSolver::sourceU(double u, double v) {
    return -u * v * v + m_f * (1 - u);
}
double GrayScottSolver::sourceV(double u, double v) {
    return u * v * v - (m_f + m_k) * v;
}


