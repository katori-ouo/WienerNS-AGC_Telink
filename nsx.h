/**
 * @file nsx.h
 * @author fjj
 * @brief 定点版本的ns头文件
 * @version 0.1
 * @date 2022-08-30
 *
 * @copyright Copyright (c) 2012 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 *
 */

#ifndef INNOTALK_NOISE_SUPRESSION_H_
#define INNOTALK_NOISE_SUPRESSION_H_

#include "audio_config.h"
#include "signal_processing_library.h"

#define CUT        0
#define LIIROPT    0
#define EngScaled  1
#define FIRAVERAGE 0
#define FSmooth    0

#include <nds_basic_math.h>
#include <nds_statistics_math.h>
#include <nds_transform_math.h>
#include <nds_utils_math.h>
#if (LIIROPT || FIRAVERAGE)
#include <nds_filtering_math.h>
#endif

#define SMOOTH_APY16	22937  // 平滑功率谱的平滑因子0.7, Q15
#define SMOOTH_APY16_S  9831   // Q15MOD - SMOOTH_APY16, Q15
#define SMOOTH_R16		32702  // 更新平滑功率谱最小值的经验常数0.998, Q15
#define SMOOTH_BETA16	31457  // 更新平滑功率谱最小值的经验常数0.96, Q15
#define AP16			6554   // 语音概率更新系数0.2, Q15
#define AP16_S          26214  // Q15MOD - AP16, Q15
#define G16				3277   // 语音概率变换系数0.1, Q15
#define G16_S           29491  // Q15MOD - G16, Q15
#define AF16			32112  // 更新先验信噪比的经验系数0.98, Q15
#define AF16_S          656    // Q15MOD - AF16, Q15
#define B_LIM16         16384  // 最终能量增益因子计算中的阈值0.5, Q15
#define FILTER_SMOOTH   29491  // 滤波器更新系数0.9, Q15
#define FILTER_SMOOTH_S 3277   // Q15MOD - FILTER_SMOOTH, Q15

#if FSmooth
#define ZETA_THR        13107   // 去除音乐噪声的信噪比阈值0.4, Q15
#define PSI             20      // 去除音乐噪声的平滑程度, Q0
#endif

#define TMPVALUE        30000   // 为了解决异常点暂时设置的值

typedef struct NsxInst_t_ {
    int                     initFlag;

    uint32_t                fs;
    int                     blockLen; // 块长度
    int16_t                 analysisBuffer[ANAL_BLOCKL_MAX];
    int16_t                 synthesisBuffer[ANAL_BLOCKL_MAX];

    int                     stages;
    int                     blockIndex;  // Frame index counter.
    int                     zeroInputSignal;  // Zero input signal flag.
    int                     scaleEnergyIn;
    uint64_t                energyIn;
    int                     norm;
    int                     normPrev;
    int16_t                 fftdata[ANAL_BLOCKL_MAX * 2];

    // 维纳滤波有关
    int16_t         smooth32[HALF_ANAL_BLOCKL];          // 滤波器系数
    int32_t         overdrive16;                         /* Q15 */
    int32_t         denoiseBound16;                      /* Q15 */
    int32_t         noisePrev32[HALF_ANAL_BLOCKL];       // 前一帧的噪声功率谱
    int32_t         pmagnPrev32[HALF_ANAL_BLOCKL];       // 前一帧的平滑功率谱
    int32_t         probPrev32[HALF_ANAL_BLOCKL];        // 前一帧的语音概率
    int32_t         signalPrev32[HALF_ANAL_BLOCKL];      // 前一帧干净语音的幅度谱估计值
    int32_t         minMagn32[HALF_ANAL_BLOCKL];         // 平滑功率谱的最小值
} NsxInst_t;
typedef struct NsxHandleT NsxHandle;

#if FSmooth
typedef struct {
    int32_t  size;
    int32_t  sum;
    int16_t* queue;
    int32_t  front;
    int32_t  rear;
} MovingAverage;
#endif

#ifdef __cplusplus
extern "C"
{
#endif
int InnoTalkNsx_Create(void** nsxInst);
int32_t InnoTalkNsx_InitCore(void* inst, uint32_t fs);
int InnoTalkNsx_ProcessCore(void* inst, short* inFrameLow, short* outFrameLow);

#if LIIROPT
int InnoTalkHpf_Create(void** hpfInst);
int32_t InnoTalkHpf_InitCore(void* inst, float* rcoeff, float* lcoeff, int nstage);
int InnoTalkHpf_ProcessFloat(void* inst, short* inFrame, short* outFrame);
#endif

#ifdef __cplusplus
}
#endif

#endif /* INNOTALK_NOISE_SUPRESSION_H_ */
