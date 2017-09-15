//#include "Backend/BackOptimize.h"
//#include "Frontend/Tracking.h"

#include "Backend/BackEnd.h"

namespace SLAMSystem {


float EnergyFunction::backGNOptimize(int IterNum)
{
    {/*    Tracking* pTracking = new Tracking();
    if(pTracking->allKeyframes.size() < 2) return 0;
    if(pTracking->allKeyframes.size() < 3) IterNum = 20;
    if(pTracking->allKeyframes.size() < 4) IterNum = 15;

    pTracking->allResiduals.clear();

    for(Frame* pFrame:pTracking->allKeyframes)
    {
        for(PixelPoint* pPoint:pFrame->pixelPointsActive)
        {
            for(Residual* pResidual:pPoint->pointResiduals)
            {
                if(!pResidual->backData->isLinearized)
                {
                    pTracking->allResiduals.push_back(pResidual);
                    pResidual->resetOOB();
                }
                else
                {
                    assert(false||"Todo");
                }
            }
        }
    }

    Vec3 lastEnergy = pTracking->linearizeAll(false);
    */}

    int numPoints = 0;
    int numLinearizedRes = 0;
    m_pTracking->allResiduals.clear();
    for(Frame* frame : m_pTracking->allKeyframes)
    {
        for(PixelPoint* point : frame->pixelPointsActive)
        {
            for(Residual* res : point->pointResiduals)
            {
                //所有关键帧中所有激活的像素点的所有的残差项如果没有线性化，那么放到allResiduals中
                if(res->backData->isLinearized == false)
                {
                    m_pTracking->allResiduals.push_back(res);
                    res->resetOOB();
                }
                else
                    numLinearizedRes++;
            }
            numPoints++;
        }
    }
    cout << "Optimize " << numPoints << " points, "
         << (int)m_pTracking->allResiduals.size() << " active residuals, "
         << numLinearizedRes << " linearized residuals." << endl;

    Vec3 lastEnergy = accumulateAllRes(false);  // 为true的情况放在后面再讨论
    double lastEnergyLinearize = calculateLinearizeEnergy();    // setting_forceAceptStep 为ｆａｌｓｅ的情况放到后边讨论
    double lastEnergyMarginalize = calculateMarginalizeEnergy();

    applyResidual();

    double lambda = 1e-1;
    float setpSize = 1;

    /* *** GN迭代优化的核心部分 *** */
    for(int iteration = 0; iteration < IterNum; iteration++)
    {
        backupState(iteration!=0);

        solveSystem(iteration, lambda);

        bool ok = addStepFromBackup();


    }
}

Vec3 EnergyFunction::accumulateAllRes(bool fixLinearization)    // todo fixLinearization
{
    double lastEnergy = 0;
    double lastEnergyR = 0; // not used
    double num = 0; // not used
    vector<Residual*> res_toRemove[ThreadNum];
    for(int i = 0; i < ThreadNum; i++)
        res_toRemove[i].clear();

    if(settings_UseMultiThread)
    {
        auto newFunc = std::bind(&EnergyFunction::accumulateAllRes_Reductor, this, fixLinearization, res_toRemove, _1, _2, _3, _4);
        threadReduce1->reduce(newFunc, 0, m_pTracking->allResiduals.size(), 0);
        lastEnergy = (double)threadReduce1->basicUnit;
    }
    else
    {
        float basicUnit;
        accumulateAllRes_Reductor(fixLinearization, res_toRemove, 0, m_pTracking->allResiduals.size(), &basicUnit, 0);
        lastEnergy = basicUnit;
    }

    /* 对于直接法，更新帧的能量阈值 */
    setNewFrameEnergyTH();
    /* 特征点法没有这一步 */

    if(fixLinearization)
    {
        // 8.23 todo    9.15 todo
    }

    return Vec3(lastEnergy, lastEnergyR, num);
}

void EnergyFunction::accumulateAllRes_Reductor(bool fixLinearization, vector<Residual *> *res_toRemove, int min, int max, float *basicUnit, int threadId)
{
    for(int i = min; i < max; i++)
    {
        Residual* res = m_pTracking->allResiduals[i];
        *basicUnit += res->calcuResidualAndJacobianPart1(m_pTracking->camera, trackMode, cameraMode);
    }

    if(fixLinearization)
    {
        // todo 8.22
    }
}

void EnergyFunction::setNewFrameEnergyTH()  // 8.23 finished
{
    // todo 8.22
    m_pTracking->allResVec.clear();
    m_pTracking->allResVec.reserve(m_pTracking->allResiduals.size()*2);

    Frame* nowFrame = m_pTracking->allKeyframes.back();
    for(Residual* r : m_pTracking->allResiduals)
    {
        if (r->newEnergyForTH >=0 && r->targetFrame == nowFrame)
            m_pTracking->allResVec.push_back(r->newEnergyForTH);
    }
    assert(m_pTracking->allResVec.size() > 0);

    int nthIdx = settings_frameEnergyTH * m_pTracking->allResVec.size();
    // 取得allResVec中元素值大小为第nthIdx大的元素值allResVec[nthIdx]
    assert(nthIdx < m_pTracking->allResVec.size());
    std::nth_element(m_pTracking->allResVec.begin(), m_pTracking->allResVec.begin()+nthIdx, m_pTracking->allResVec.end());
    float sqrtNthElement = sqrt(m_pTracking->allResVec[nthIdx]);

    // 设置frameEnergyTH
    nowFrame->frameEnergyTH = sqrtNthElement*settings_frameEnergyTHFacMedian;
    nowFrame->frameEnergyTH = 26.0f*settings_frameEnergyTHConstWeight + nowFrame->frameEnergyTH*(1-settings_frameEnergyTHConstWeight);
    nowFrame->frameEnergyTH = nowFrame->frameEnergyTH * nowFrame->frameEnergyTH;
    nowFrame->frameEnergyTH *= settings_overallEnergyTHWeight*settings_overallEnergyTHWeight;

}

double EnergyFunction::calculateLinearizeEnergy()
{
    if(settings_forceAceptStep) return 0;
    else
    {
        // false todo 8.23
    }
}

double EnergyFunction::calculateMarginalizeEnergy()
{
    if(settings_forceAceptStep) return 0;
    {
        // false todo 8.23
    }
}

void EnergyFunction::applyResidual()    // 8.23 finished
{
    if(settings_UseMultiThread)
    {
        auto newFunc = std::bind(&EnergyFunction::applyResidual_Reductor, this, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc, 0, m_pTracking->allResiduals.size(), 50);
    }
    else
        applyResidual_Reductor(0, m_pTracking->allResiduals.size(), NULL, 0);
}

void EnergyFunction::applyResidual_Reductor(int min, int max, Vec10 *basicUnit, int threadId)
{
    for(int i = min; i < max; i++)
    {
        Residual* res = m_pTracking->allResiduals[i];
        res->applyRes();
    }
}






void EnergyFunction::backupState(bool backupLastStep)   // finished 8.24
{
    if(backupLastStep)
    {
        m_pTracking->camera->backData->step_backup = m_pTracking->camera->backData->step;
        m_pTracking->camera->backData->intriParas_backup = m_pTracking->camera->backData->intriParas;
        for(Frame* frame : m_pTracking->allKeyframes)
        {
            frame->backData->step_backup = frame->backData->step;
            frame->backData->state_x_backup = frame->backData->state_x;
            for(PixelPoint* point : frame->pixelPointsActive)
            {
                point->backData->idepth_backup = point->backData->idepth;
                point->backData->step_backup = point->backData->step;
            }
        }
    }
    else    // 第一次优化
    {
        m_pTracking->camera->backData->step_backup = Vec4f::Constant(0);
        m_pTracking->camera->backData->intriParas_backup = m_pTracking->camera->backData->intriParas;
        for(Frame* frame : m_pTracking->allKeyframes)
        {
/* 特征点法和直接法的区别 Vec8f Vec6f */
            frame->backData->step_backup = Vec8::Constant(0);
            frame->backData->state_x_backup = frame->backData->state_x;
            for(PixelPoint* point : frame->pixelPointsActive)
            {
                point->backData->idepth_backup = point->backData->idepth;
                point->backData->step_backup = 0;
            }
        }
    }
}

void EnergyFunction::solveSystem(int iteration, double lambda)
{
    getNullspace(lastNullspace_pose, lastNullspace_scale, lastNullspace_affA, lastNullspace_affB);
    this->solve(iteration, lambda, m_pTracking->camera);
}


void EnergyFunction::getNullspace(vector<VecX> &nullspace_pose, vector<VecX> &nullspace_scale, vector<VecX> &nullspace_affA, vector<VecX> &nullspace_affB)
{
    nullspace_pose.clear();
    nullspace_scale.clear();
    nullspace_affA.clear(); /* 特征点法没有光度参数 */
    nullspace_affB.clear();

    /* 直接法，且优化相机内参时的n */
    {
        int n = CPARS + m_pTracking->allKeyframes.size() * 8;
        for(int i = 0; i < 6; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();

            for(Frame* frame : m_pTracking->allKeyframes)
            {
                nullspace_x0.segment<6>(CPARS + frame->keyframeId * 8) = frame->backData->nullspace_pose.col(i);
            }
            nullspace_pose.push_back(nullspace_x0);
        }

        for(int i = 0; i < 2; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();
            for(Frame* frame : m_pTracking->allKeyframes)
            {
                nullspace_x0.segment<2>(CPARS + frame->keyframeId * 8 + 6) = frame->backData->nullspace_affine.col(i).head<2>();    // (a,b)
            }
            if(i==0) nullspace_affA.push_back(nullspace_x0);
            if(i==1) nullspace_affB.push_back(nullspace_x0);
        }

        VecX nullspace_x0(n); nullspace_x0.setZero();
        for(Frame* frame : m_pTracking->allKeyframes)
        {
            nullspace_x0.segment<6>(CPARS + frame->keyframeId * 8) = frame->backData->nullspace_scale;
        }
        nullspace_scale.push_back(nullspace_x0);
    }

    /* 直接法，不优化相机内参时的n */
    {
        int n = m_pTracking->allKeyframes.size() * 8;

        for(int i = 0; i < 6; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();

            for(Frame* frame : m_pTracking->allKeyframes)
            {
                nullspace_x0.segment<6>(frame->keyframeId * 8) = frame->backData->nullspace_pose.col(i);
            }
            nullspace_pose.push_back(nullspace_x0);
        }

        for(int i = 0; i < 2; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();
            for(Frame* frame : m_pTracking->allKeyframes)
            {
                nullspace_x0.segment<2>(frame->keyframeId * 8 + 6) = frame->backData->nullspace_affine.col(i).head<2>();    // (a,b)
            }
            if(i==0) nullspace_affA.push_back(nullspace_x0);
            if(i==1) nullspace_affB.push_back(nullspace_x0);
        }

        VecX nullspace_x0(n); nullspace_x0.setZero();
        for(Frame* frame : m_pTracking->allKeyframes)
        {
            nullspace_x0.segment<6>(frame->keyframeId * 8) = frame->backData->nullspace_scale;
        }
        nullspace_scale.push_back(nullspace_x0);
    }
    /* 特征点法，优化相机内参时的n */
    {
        int n = CPARS + m_pTracking->allKeyframes.size() * 6;

        for(int i = 0; i < 6; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();

            for(Frame* frame : m_pTracking->allKeyframes)
            {
                nullspace_x0.segment<6>(CPARS + frame->keyframeId * 6) = frame->backData->nullspace_pose.col(i);
            }
            nullspace_pose.push_back(nullspace_x0);
        }

        // 特征点法(a,b)为0
        for(int i = 0; i < 2; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();
            // frame->backData->nullspace_affine.col(i).head<2>()   特征点法的nullspace_affine全部设置为０
            if(i==0) nullspace_affA.push_back(nullspace_x0);
            if(i==1) nullspace_affB.push_back(nullspace_x0);
        }

        // 特征点法貌似不需要scale之类的东西
        VecX nullspace_x0(n); nullspace_x0.setZero();
        for(Frame* frame : m_pTracking->allKeyframes)
        {
            nullspace_x0.segment<6>(CPARS + frame->keyframeId * 6) = frame->backData->nullspace_scale;  //frame->backData->nullspace_scale  特征点法的nullspace_scale全部设置为１
        }
        nullspace_scale.push_back(nullspace_x0);
    }
    /* 特征点法，不优化相机内参时的n */
    {
        int n = m_pTracking->allKeyframes.size() * 6;
        for(int i = 0; i < 6; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();

            for(Frame* frame : m_pTracking->allKeyframes)
            {
                nullspace_x0.segment<6>(frame->keyframeId * 6) = frame->backData->nullspace_pose.col(i);
            }
            nullspace_pose.push_back(nullspace_x0);
        }

        for(int i = 0; i < 2; i++)
        {
            VecX nullspace_x0(n); nullspace_x0.setZero();
            // frame->backData->nullspace_affine.col(i).head<2>()   特征点法的nullspace_affine全部设置为０
            if(i==0) nullspace_affA.push_back(nullspace_x0);
            if(i==1) nullspace_affB.push_back(nullspace_x0);
        }

        VecX nullspace_x0(n); nullspace_x0.setZero();
        for(Frame* frame : m_pTracking->allKeyframes)
        {
            nullspace_x0.segment<6>(frame->keyframeId * 6) = frame->backData->nullspace_scale;
            //frame->backData->nullspace_scale  特征点法的nullspace_scale全部设置为１
        }
        nullspace_scale.push_back(nullspace_x0);

    }
}

/**
 * @brief EnergyFunction::solve         得到H,b,求解增量deltaX，然后将增量返回到相应的状态量中
 * @param iteration
 * @param lambda
 * @param cam
 */
void EnergyFunction::solve(int iteration, double lambda, Camera *cam)   // todo 8.25 afternoon
{
    if(settings_solverMode & SOLVER_USE_GN) lambda=0;
    if(settings_solverMode & SOLVER_FIX_LAMBDA) lambda = 1e-5;

    MatXX H_activate_top, H_linear_top, H_schur;
    VecX  b_activate_top, b_linear_top, b_schur;
    VecX b_marg_top;

    accumulateActivatePart(H_activate_top, b_activate_top, settings_UseMultiThread);

    accumulateLinearPart(H_linear_top, b_linear_top, settings_UseMultiThread);

    accumulateSchurPart(H_schur, b_schur, settings_UseMultiThread);

    b_marg_top = (b_M + H_M * getStitchDelta());    // 先求出帧的边缘化部分的b

    MatXX H_Final_top;
    VecX b_Final_top;

    H_Final_top = H_linear_top + H_activate_top + H_M - H_schur;
    b_Final_top = b_linear_top + b_activate_top + b_marg_top - b_schur;     // ??? 为什么不是b_M

    last_H = H_Final_top;
    last_b = b_Final_top;

    for(int i = 0; i < H_Final_top.cols(); i++)
        H_Final_top(i, i)  *= (1 + lambda);

    calculateDeltaXBySVD(last_H, last_b, last_deltaX);

    feedBackDeltaX(last_deltaX, cam, settings_UseMultiThread);

}


void EnergyFunction::accumulateActivatePart(MatXX &H, VecX &b, bool MT)
{
    if(settings_UseMultiThread)
    {
        // 初始化
        auto newFunc = std::bind(&AccumulatedTopHessian::setZero, acc_top_Activate, allBackKeyFramesNum, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc, 0, allBackKeyFramesNum, 0);
        // 计算刚激活的点的残差及H,b
        auto newFunc2 = std::bind(&AccumulatedTopHessian::addPointMultiThread<active>, acc_top_Activate, &allBackPixelPoints, this, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc2, 0, allBackPixelPointsNum, 0);
        // 得到全局的H,b
        acc_top_Activate->stitchDoubleMultiThread(threadReduce2, H, b, this, false, true);
        resNum_Active = acc_top_Activate->resNum[0];
    }
    else
    {
        // todo 9.7
    }
}

void EnergyFunction::accumulateLinearPart(MatXX &H, VecX &b, bool MT)
{
    if(settings_UseMultiThread)
    {
        // 初始化
        auto newFunc = std::bind(&AccumulatedTopHessian::setZero, acc_top_Linearized, allBackKeyFramesNum, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc, 0, allBackKeyFramesNum, 0);
        // 计算线性化的点的残差及H,b
        auto newFunc2 = std::bind(&AccumulatedTopHessian::addPointMultiThread<linearized>, acc_top_Linearized, &allBackPixelPoints, this, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc2, 0, allBackKeyFramesNum, 0);
        // 得到全局的H,b
        acc_top_Linearized->stitchDoubleMultiThread(threadReduce2, H, b, this, true, true);
        resNum_Linearized = acc_top_Linearized->resNum[0];
    }
    else
    {
        // todo 9.7
    }
}

void EnergyFunction::accumulateSchurPart(MatXX &H, VecX &b, bool MT)
{

    if(settings_UseMultiThread)
    {
        // 初始化
        auto newFunc = std::bind(&AccumulatedSCHessian::setZero, acc_bot_Marginalize, allBackKeyFramesNum, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc, 0, allBackKeyFramesNum, 0);
        // 计算边缘化的点的残差及该部分的H,b
        auto newFunc2 = std::bind(&AccumulatedSCHessian::addPointMultiThread, acc_bot_Marginalize, &allBackPixelPoints, true, _1, _2, _3, _4);
        threadReduce2->reduce(newFunc2, 0, allBackPixelPointsNum, 0);
        // 得到全局的H,b
        acc_bot_Marginalize->stitchDoubleMultiThread(threadReduce2, H, b, this, true);
    } else
    {
        acc_bot_Marginalize->setZero(allBackKeyFramesNum);
        for(BackFrame* frame : allBackKeyFrames)
        {
            for(BackPixelPoint* point : frame->backPixelPoints)
            {
                acc_bot_Marginalize->addPoint(point, true);
            }
        }
        acc_bot_Marginalize->stitchDoubleMultiThread(threadReduce2, H, b, this, false);
    }
}

/**
 * @brief 得到增量delta_top(不包含idepth的增量)
 *
 */
VecX EnergyFunction::getStitchDelta()
{
    if(trackMode == Direct && cameraMode == Optimize)
    {
        VecX delta_top = VecX(CPARS + allBackKeyFramesNum * 8);
        delta_top.head<CPARS>() = deltaF_C.cast<double>();

        for(int i = 0; i < allBackKeyFramesNum; i++)
        {
            delta_top.segment<8>(CPARS + 8 * i) = allBackKeyFrames[i]->delta_0ToNow;
        }
        return delta_top;
    } else if (trackMode == Direct && cameraMode == Fix)
    {
        assert(deltaF_C == Vec4f::Zero());

        VecX delta_top = VecX(allBackKeyFramesNum * 8);
        for(int i = 0; i < allBackKeyFramesNum; i++)
        {
            delta_top.segment<8>(8 * i) = allBackKeyFrames[i]->delta_0ToNow;
        }
        return delta_top;
    } else if(trackMode == Feature && cameraMode == Optimize)
    {
        VecX delta_top = VecX(CPARS + allBackKeyFramesNum * 6);
        delta_top.head<CPARS>() = deltaF_C.cast<double>();

        for(int i = 0; i < allBackKeyFramesNum; i++)
        {
            delta_top.segment<6>(CPARS + 6 * i) = allBackKeyFrames[i]->delta_0ToNow_Feature;
        }
        return delta_top;
    } else if(trackMode == Feature && cameraMode == Fix)
    {
        assert(deltaF_C == Vec4f::Zero());

        VecX delta_top = VecX(allBackKeyFramesNum * 6);
        for(int i = 0; i < allBackKeyFramesNum; i++)
        {
            delta_top.segment<6>(6 * i) = allBackKeyFrames[i]->delta_0ToNow_Feature;
        }
        return delta_top;
    } else { assert(false && "No such mode!"); }
}

/**
 * @brief 纯数学计算模块, 利用SVD分解求解线性增量方程的deltaX
 * @param deltaX
 */
void EnergyFunction::calculateDeltaXBySVD(MatXX &H, VecX &b, VecX &deltaX)
{
    VecX SVecI = H.diagonal().cwiseSqrt().cwiseInverse();
    MatXX H_Scaled = SVecI.asDiagonal() * H * SVecI.asDiagonal();
    VecX b_Scaled = SVecI.asDiagonal() * b;
    Eigen::JacobiSVD<MatXX> svd(H_Scaled, Eigen::ComputeThinU | Eigen::ComputeThinV);

    VecX S = svd.singularValues();
    double minSv = 1e10, maxSv = 0;
    for(int i = 0; i < S.size(); i++)
    {
        if(S[i] < minSv) minSv = S[i];
        if(S[i] > maxSv) maxSv = S[i];
    }

    VecX Ub = svd.matrixU().transpose() * b_Scaled;
    int setZero = 0;
    for(int i = 0; i < Ub.size(); i++)
    {
        Ub[i] /= S[i];
    }


    // 得到deltaｘ
    deltaX = SVecI.asDiagonal() * svd.matrixV() * Ub;
}

/**
 * @brief EnergyFunction::feedBackDeltaX        deltaｘ作用到c,ab,xi,的step中,得到各个状态量每一步的增量
 * @param delta_X
 * @param cam
 * @param MT
 */
void EnergyFunction::feedBackDeltaX(VecX delta_X, Camera *cam, bool MT)     // finished 09.14 18:09
{
    if(trackMode == Direct && cameraMode == Optimize)
    {
        assert(delta_X.size() == CPARS + allBackKeyFramesNum * 8);
// 1. 更新Camera的step
        VecXf delta_X_float = delta_X.cast<float>();
        cam->backData->step = - delta_X_float.head<CPARS>(); // deltaX作用到c的step中
// 2. 更新Frame的step,计算出伴随矩阵的step用于第三步的计算
        Mat18f* delta_X_xi_adjoint = new Mat18f [allBackKeyFramesNum * allBackKeyFramesNum];
        Vec4f delta_X_cam = delta_X_float.head<CPARS>();       // 用于更新点的深度的变量
        Mat88f* adHost_float = new Mat88f [allBackKeyFramesNum * allBackKeyFramesNum];
        Mat88f* adTarget_float = new Mat88f [allBackKeyFramesNum * allBackKeyFramesNum];

        for(BackFrame* hostframe : allBackKeyFrames)
        {
            hostframe->step = - delta_X.segment<8>(CPARS + 8 * hostframe->keyframeId); // deltaX作用到xi,(a,b)的step中

            for(BackFrame* targetframe : allBackKeyFrames)
            {
                adHost_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adHost[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();
                adTarget_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adTarget[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();

                delta_X_xi_adjoint[allBackKeyFramesNum*hostframe->keyframeId + targetframe->keyframeId] =
                        delta_X_float.segment<8>(CPARS+8*hostframe->keyframeId).transpose() *   adHost_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId]
                        + delta_X_float.segment<8>(CPARS+8*targetframe->keyframeId).transpose() * adTarget_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId]; // 求解delta_X_xi_adjoint, delta_X_xi_adjoint用于求解点的逆深度的step
            }
        }
// 3.  更新Point的step
        if(MT)
        {
            auto newFunc = std::bind(&EnergyFunction::feedBackDeltaIdepthDirect, this, delta_X_cam, delta_X_xi_adjoint, _1, _2, _3, _4);
            threadReduce2->reduce(newFunc, 0, allBackPixelPoints.size(), 50);
        } else
        {
            feedBackDeltaIdepthDirect(delta_X_cam, delta_X_xi_adjoint, 0, allBackPixelPoints.size(), NULL, 0);
        }
        delete [] delta_X_xi_adjoint;
        delete [] adHost_float;
        delete [] adTarget_float;

    } else if(trackMode == Direct && cameraMode == Fix)
    {
        assert(delta_X.size() == allBackKeyFramesNum * 8);
// 1. 更新Camera的step
        VecXf delta_X_float = delta_X.cast<float>();
        cam->backData->step = Vec4f::Zero();
// 2. 更新Frame的step
        Mat18f* delta_X_xi_adjoint = new Mat18f [allBackKeyFramesNum * allBackKeyFramesNum];
        Mat88f* adHost_float = new Mat88f [allBackKeyFramesNum * allBackKeyFramesNum];
        Mat88f* adTarget_float = new Mat88f [allBackKeyFramesNum * allBackKeyFramesNum];

        for(BackFrame* hostframe : allBackKeyFrames)
        {
            hostframe->step = - delta_X.segment<8>(8 * hostframe->keyframeId);    // deltaX作用到xi,(a,b)的step中

            for(BackFrame* targetframe : allBackKeyFrames)
            {
                adHost_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adHost[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();
                adTarget_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adTarget[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();

                delta_X_xi_adjoint[allBackKeyFramesNum*hostframe->keyframeId + targetframe->keyframeId] =
                        delta_X_float.segment<8>(8*hostframe->keyframeId).transpose() * adHost_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId]
                        + delta_X_float.segment<8>(8*targetframe->keyframeId).transpose() * adTarget_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId]; // 求解delta_X_xi_adjoint, delta_X_xi_adjoint用于求解点的逆深度的step
            }
        }
// 3. 更新Point的step
        if(MT)
        {
            auto newFunc = std::bind(&EnergyFunction::feedBackDeltaIdepthDirect, this, Vec4f::Zero(), delta_X_xi_adjoint, _1, _2, _3, _4);
            threadReduce2->reduce(newFunc, 0, allBackPixelPoints.size(), 50);
        } else
        {
            feedBackDeltaIdepthDirect(Vec4f::Zero(), delta_X_xi_adjoint, 0, allBackPixelPoints.size(), NULL, 0);
        }
        delete [] delta_X_xi_adjoint;
        delete [] adHost_float;
        delete [] adTarget_float;
    } else if(trackMode == Feature && cameraMode == Optimize)
    {
        assert(delta_X.size() == CPARS + allBackKeyFramesNum * 6);  // 确定输入deltaX的维度是否正确
// 1. 更新Camera的step
        VecXf delta_X_float = delta_X.cast<float>();
        cam->backData->step = - delta_X_float.head<CPARS>(); // deltaX作用到c的step中
// 2. 更新Frame的step,计算出伴随矩阵的step用于第三步的计算
        Mat16f* delta_X_xi_adjoint = new Mat16f [allBackKeyFramesNum * allBackKeyFramesNum];
        Vec4f delta_X_cam = delta_X_float.head<CPARS>();    // 用于更新点的深度的变量
        Mat66f* adHost_Feature_float = new Mat66f [allBackKeyFramesNum * allBackKeyFramesNum];
        Mat66f* adTarget_Feature_float = new Mat66f [allBackKeyFramesNum * allBackKeyFramesNum];

        for(BackFrame* hostframe : allBackKeyFrames)
        {
            hostframe->step_Feature = - delta_X.segment<6>(CPARS + 6 * hostframe->keyframeId);// deltaX作用到xi的step中, 没有光度参数(a,b)

            for(BackFrame* targetframe : allBackKeyFrames)
            {
                adHost_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adHost_Feature[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();
                adTarget_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adTarget_Feature[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();

                delta_X_xi_adjoint[allBackKeyFramesNum*hostframe->keyframeId + targetframe->keyframeId] =
                        delta_X_float.segment<6>(CPARS+6*hostframe->keyframeId).transpose() * adHost_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId]
                        + delta_X_float.segment<6>(CPARS+6*targetframe->keyframeId).transpose() * adTarget_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId];
            }
        }
// 3. 更新Point的step
        if(MT)
        {
            auto newFunc = std::bind(&EnergyFunction::feedBackDeltaIdepthFeature, this, delta_X_cam, delta_X_xi_adjoint, _1, _2, _3, _4);
            threadReduce2->reduce(newFunc, 0, allBackPixelPoints.size(), 50);
        } else
        {
            feedBackDeltaIdepthFeature(delta_X_cam, delta_X_xi_adjoint, 0, allBackPixelPoints.size(), NULL, 0);
        }
        delete [] delta_X_xi_adjoint;
        delete [] adHost_Feature_float;
        delete [] adTarget_Feature_float;
    } else if(trackMode == Feature && cameraMode == Fix)
    {
        assert(delta_X.size() == allBackKeyFramesNum * 6);
// 1. 更新Camera的step
        VecXf delta_X_float = delta_X.cast<float>();
        cam->backData->step = Vec4f::Zero();
// 2. 更新Frame的step
        Mat16f* delta_X_xi_adjoint = new Mat16f [allBackKeyFramesNum * allBackKeyFramesNum];
        Mat66f* adHost_Feature_float = new Mat66f [allBackKeyFramesNum * allBackKeyFramesNum];
        Mat66f* adTarget_Feature_float = new Mat66f [allBackKeyFramesNum * allBackKeyFramesNum];

        for(BackFrame* hostframe : allBackKeyFrames)
        {
            hostframe->step_Feature = - delta_X.segment<6>(6 * hostframe->keyframeId);

            for(BackFrame* targetframe : allBackKeyFrames)
            {
                adHost_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adHost_Feature[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();
                adTarget_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId] = adTarget_Feature[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId].cast<float>();

                delta_X_xi_adjoint[allBackKeyFramesNum*hostframe->keyframeId + targetframe->keyframeId] =
                        delta_X_float.segment<6>(6*hostframe->keyframeId).transpose() * adHost_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId]
                        + delta_X_float.segment<6>(6*targetframe->keyframeId).transpose() * adTarget_Feature_float[hostframe->keyframeId+allBackKeyFramesNum*targetframe->keyframeId];
            }
        }
// 3. 更新Point的step
        if(MT)
        {
            auto newFunc = std::bind(&EnergyFunction::feedBackDeltaIdepthFeature, this, Vec4f::Zero(), delta_X_xi_adjoint, _1, _2, _3, _4);
            threadReduce2->reduce(newFunc, 0, allBackPixelPoints.size(), 50);
        } else
        {
            feedBackDeltaIdepthFeature(Vec4f::Zero(), delta_X_xi_adjoint, 0, allBackPixelPoints.size(), NULL, 0);
        }
        delete [] delta_X_xi_adjoint;
        delete [] adHost_Feature_float;
        delete [] adTarget_Feature_float;
    } else { assert(false && "No such mode!"); }
}

void EnergyFunction::feedBackDeltaIdepthDirect(const Vec4f &deltaX_cam, Mat18f *deltaX_xi_adjoint, int min, int max, Vec10 *basicUnit, int thd_idx)
{
    assert(trackMode == Direct);

    if(cameraMode == Optimize)       // to be checked   09.14 afternoon
    {
        for(int i = min; i < max; i++)
        {
            BackPixelPoint* point = allBackPixelPoints[i];
            int goodresidual = 0;
            for(BackResidual* res : point->backResiduals)   if(res->isActiveAndIsGoodNEW) goodresidual++;
            if(goodresidual == 0)
            {
                point->step = 0;
                continue;
            }

            //******* b = bdSum - xc*Hc - xigema(xAd(hi,ti) * JpJd) (1.188)
            float b = point->b_Sum - deltaX_cam.dot(point->H_C_idepth_accAF + point->H_C_idepth_accLF);
            for(BackResidual* res : point->backResiduals)
            {
                if(res->isActiveAndIsGoodNEW == false) continue;
                b = b - deltaX_xi_adjoint[res->backHostIndex * allBackKeyFramesNum + res->backTargetIndex] * res->JpJdF;
            }
            point->step = - b * point->H_depth;
            assert(std::isfinite(point->step));
        }
    } else if(cameraMode == Fix)
    {
        assert(deltaX_cam == Vec4f::Zero());

        for(int i = min; i < max; i++)
        {
            BackPixelPoint* point = allBackPixelPoints[i];
            int goodresidual = 0;
            for(BackResidual* res : point->backResiduals)   if(res->isActiveAndIsGoodNEW) goodresidual++;
            if(goodresidual == 0)
            {
                point->step = 0;
                continue;
            }

            //******* b = bdSum - xigma(delta_adjoint(hi,ti) * JpJd) (1.188)    由于deltaX_cam为0
            float b = point->b_Sum;
            for(BackResidual* res : point->backResiduals)
            {
                if(res->isActiveAndIsGoodNEW == false) continue;
                b = b - deltaX_xi_adjoint[res->backHostIndex * allBackKeyFramesNum + res->backTargetIndex] * res->JpJdF;
            }
            point->step = - b * point->H_depth;
            assert(std::isfinite(point->step));
        }
    } else { assert(false && "No such mode!"); }
}

void EnergyFunction::feedBackDeltaIdepthFeature(const Vec4f &deltaX_cam, Mat16f *deltaX_xi_adjoint, int min, int max, Vec10 *basicUnit, int thd_idx)
{
    assert(trackMode == Feature);

    if(cameraMode == Optimize)
    {
        for(int i = min; i < max; i++)
        {
            BackPixelPoint* point = allBackPixelPoints[i];
            int goodresidual = 0;
            for(BackResidual* res : point->backResiduals)   if(res->isActiveAndIsGoodNEW)   goodresidual++;
            if(goodresidual == 0)
            {
                point->step = 0;
                continue;
            }

            //******** b = bdSum - deltaX_cam*Hc - xigma(delta_adjoint(hi,ti) * JpJd) (1.188)
            float b = point->b_Sum - deltaX_cam.dot(point->H_C_idepth_accAF + point->H_C_idepth_accLF);
            for(BackResidual* res : point->backResiduals)
            {
                if(res->isActiveAndIsGoodNEW == false) continue;
                b = b - deltaX_xi_adjoint[res->backHostIndex * allBackKeyFramesNum + res->backTargetIndex] * res->JpJdF_Feature;
            }
            point->step = - b * point->H_depth;
            assert(std::isfinite(point->step));
        }
    } else if(cameraMode == Fix)
    {
        for(int i = min; i < max; i++)
        {
            BackPixelPoint* point = allBackPixelPoints[i];
            int goodresidual = 0;
            for(BackResidual* res : point->backResiduals)   if(res->isActiveAndIsGoodNEW)   goodresidual++;
            if(goodresidual == 0)
            {
                point->step = 0;
                continue;
            }

            float b = point->b_Sum - deltaX_cam.dot(point->H_C_idepth_accAF + point->H_C_idepth_accLF);
            for(BackResidual* res : point->backResiduals)
            {
                if(res->isActiveAndIsGoodNEW == false) continue;
                b = b - deltaX_xi_adjoint[res->backHostIndex * allBackKeyFramesNum + res->backTargetIndex] * res->JpJdF_Feature;
            }
            point->step = - b * point->H_depth;
            assert(std::isfinite(point->step));
        }
    } else { assert(false && "No such mode!"); }
}

bool EnergyFunction::addStepFromBackup()
{
    float sumT = 0, sumR = 0;   // R,T
    float sumA = 0, sumB = 0;   // (a,b)
    float sumPnt = 0, sumNID = 0;   // point_step, point_idepth
    float numPoint = 0;
// 1. 分情况更新状态量中的相机内参
    if(cameraMode == Optimize)
    {
        // x_new = x + deltax
        m_pTracking->camera->backData->intriParas = m_pTracking->camera->backData->intriParas_backup + m_pTracking->camera->backData->step;
    } else if(cameraMode == Fix)
    {
        assert(m_pTracking->camera->backData->intriParas == m_pTracking->camera->backData->intriParas_0);
    }

    if(trackMode == Direct)
    {
        for(BackFrame* frame : allBackKeyFrames)
        {
// 2. 帧的状态量xi,(a,b)相加，对应的李代数位姿xi相乘
            frame->state_x = frame->state_x_backup + frame->step;
            frame->Tcw = SE3::exp(frame->state_x.segment<6>(0)) * frame->Tcw0;
            frame->Twc = frame->Tcw.inverse();

            sumT += frame->step.segment<3>(0).squaredNorm();
            sumR += frame->step.segment<3>(3).squaredNorm();
            sumA += frame->step[6] * frame->step[6];
            sumB += frame->step[7] * frame->step[7];

// 3. 点的状态量逆深度(idepth)的更新
            for(BackPixelPoint* point : frame->backPixelPoints)
            {
                point->idepth = point->idepth_backup + point->step;
                // 为什么要改变点的逆深度初始值？？？ --> 为了每次的point->delta的计算
                point->idepth0 = point->idepth_backup + point->step;

                sumPnt += point->step * point->step;
                sumNID += fabs(point->idepth_backup);
                numPoint++;
            }
        }

        m_pTracking->setPreCalcValues();

        sumA /= allBackKeyFrames.size();
        sumB /= allBackKeyFrames.size();
        sumR /= allBackKeyFrames.size();
        sumT /= allBackKeyFrames.size();
        sumPnt /= numPoint;
        sumNID /= numPoint;

        return (sqrt(sumT) < 0.0005 * settings_thresGNOptimize &&
                sqrt(sumR) < 0.0005 * settings_thresGNOptimize &&
                sqrt(sumA) < 0.0005 * settings_thresGNOptimize &&
                sqrt(sumB) < 0.0005 * settings_thresGNOptimize);

    }
    else if(trackMode == Feature)
    {
        for(BackFrame* frame : allBackKeyFrames)
        {
// 2. 帧的状态量xi,(a,b)相加，对应的李代数位姿xi相乘
            frame->state_x_Feature = frame->state_x_backup_Feature + frame->step_Feature;
            frame->Tcw = SE3::exp(frame->state_x_Feature) * frame->Tcw0;
            frame->Twc = frame->Tcw.inverse();

            sumT += frame->step.segment<3>(0).squaredNorm();
            sumR += frame->step.segment<3>(3).squaredNorm();

// 3. 点的状态量逆深度(idepth)的更新
            for(BackPixelPoint* point : frame->backPixelPoints)
            {
                point->idepth = point->idepth_backup + point->step;
                // 为什么要改变点的逆深度初始值？？？ --> 为了每次的point->delta的计算
                point->idepth0 = point->idepth_backup + point->step;

                sumPnt += point->step * point->step;
                sumNID += fabs(point->idepth_backup);
                numPoint++;
            }
        }

        m_pTracking->setPreCalcValues();

        sumA /= allBackKeyFrames.size();
        sumB /= allBackKeyFrames.size();
        sumR /= allBackKeyFrames.size();
        sumT /= allBackKeyFrames.size();
        sumPnt /= numPoint;
        sumNID /= numPoint;

        return (sqrt(sumT) < 0.0005 * settings_thresGNOptimize &&
                sqrt(sumR) < 0.0005 * settings_thresGNOptimize &&
                sqrt(sumA) < 0.0005 * settings_thresGNOptimize &&
                sqrt(sumB) < 0.0005 * settings_thresGNOptimize);

    }
}


}
