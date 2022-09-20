/**
 * @file nsx.h
 * @author fjj
 * @brief ����汾��nsͷ�ļ�
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

#ifndef NSX_H_
#define NSX_H_

#include "signal_processing_library.h"

#define  Q15MOD               32768
#define  Q15MAX               32767
#define  Q31MOD          2147483648
#define  Q30MOD          1073741824

#define CEVAOPT    0
#define TelinkOPT  1
#define CUT        0
#define LIIROPT    0
#define EngScaled  1
#define FIRAVERAGE 0
#define FSmooth    0

#if CEVAOPT
#include "./cevalib/CEVA_common.h"
#include "./cevalib/ceva_typedef.h"
#endif

#if TelinkOPT
#include <nds_basic_math.h>
#include <nds_statistics_math.h>
#include <nds_transform_math.h>
#include <nds_utils_math.h>
#if (LIIROPT || FIRAVERAGE)
#include <nds_filtering_math.h>
#endif
#endif

#define FRAME_LEN               128 /* frame length */
#define ANAL_BLOCKL_MAX         256 /* Max analysis block length */
#define ANAL_BLOCKL_MAXD2       128 /* Max analysis block length */
#define HALF_ANAL_BLOCKL        129 /* Half max analysis block length + 1 */
#define END_STARTUP_LONG        200

#define SMOOTH_APY16	22937  // ƽ�������׵�ƽ������0.7, Q15
#define SMOOTH_APY16_S  9831   // Q15MOD - SMOOTH_APY16, Q15
#define SMOOTH_R16		32702  // ����ƽ����������Сֵ�ľ��鳣��0.998, Q15
#define SMOOTH_BETA16	31457  // ����ƽ����������Сֵ�ľ��鳣��0.96, Q15
#define AP16			6554   // �������ʸ���ϵ��0.2, Q15
#define AP16_S          26214  // Q15MOD - AP16, Q15
#define G16				3277   // �������ʱ任ϵ��0.1, Q15
#define G16_S           29491  // Q15MOD - G16, Q15
#define AF16			32112  // ������������ȵľ���ϵ��0.98, Q15
#define AF16_S          656    // Q15MOD - AF16, Q15
#define B_LIM16         16384  // ���������������Ӽ����е���ֵ0.5, Q15
#define FILTER_SMOOTH   29491  // �˲�������ϵ��0.9, Q15
#define FILTER_SMOOTH_S 3277   // Q15MOD - FILTER_SMOOTH, Q15

#if FSmooth
#define ZETA_THR        13107   // ȥ�������������������ֵ0.4, Q15
#define PSI             20      // ȥ������������ƽ���̶�, Q0
#endif

#define TMPVALUE        30000   // Ϊ�˽���쳣����ʱ���õ�ֵ

typedef struct NsxInst_t_ {
    int                     initFlag;

    uint32_t                fs;
    int                     blockLen; // �鳤��
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

    // ά���˲��й�
    int16_t         smooth32[HALF_ANAL_BLOCKL];          // �˲���ϵ��
    int32_t         overdrive16;                         /* Q15 */
    int32_t         denoiseBound16;                      /* Q15 */
    int32_t         noisePrev32[HALF_ANAL_BLOCKL];       // ǰһ֡������������
    int32_t         pmagnPrev32[HALF_ANAL_BLOCKL];       // ǰһ֡��ƽ��������
    int32_t         probPrev32[HALF_ANAL_BLOCKL];        // ǰһ֡����������
    int32_t         signalPrev32[HALF_ANAL_BLOCKL];      // ǰһ֡�ɾ������ķ����׹���ֵ
    int32_t         minMagn32[HALF_ANAL_BLOCKL];         // ƽ�������׵���Сֵ
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

#endif /* NSX_H_ */
