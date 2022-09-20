/**
 * @file nsx.c
 * @author fjj
 * @brief
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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nsx.h"

int16_t FFTIOtT[ANAL_BLOCKL_MAX * 2];

const int32_t delta32[129] = {
  42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598,
  42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598,
  42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598, 42598,
  42598, 42598, 41871, 43765, 45660, 47554, 49448, 51343, 53237, 55132,
  57026, 58920, 60815, 62709, 64604, 66498, 68392, 70287, 72181, 74076,
  75970, 77864, 79759, 81653, 83548, 85442, 87336, 89231, 91125, 93020,
  94914, 96808, 98703, 100597, 102492, 104386, 106280, 108175, 110069, 111964,
  113858, 115752, 117647, 119541, 121436, 123330, 125224, 127119, 129013, 130908,
  132802, 134696, 136591, 138485, 140380, 142274, 144168, 146063, 147957, 149852,
  151746, 153640, 155535, 157429, 159324, 161218, 163840, 163840, 163840, 163840,
  163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840,
  163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840,
  163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840};

static int16_t kBlocks128w256x[256] = {
  0,404,802,1206,1608,2008,2412,2812,3212,3612,4010,4410,4808,5206,5604,5996,6394,
  6786,7180,7572,7962,8352,8740,9126,9512,9896,10280,10660,11040,11416,11794,12166,
  12540,12910,13278,13644,14012,14372,14732,15090,15446,15800,16152,16498,16846,17190,
  17530,17868,18206,18536,18868,19196,19520,19842,20158,20476,20788,21096,21404,21706,
  22006,22302,22594,22886,23170,23452,23730,24010,24282,24546,24812,25074,25330,25582,
  25832,26076,26320,26558,26792,27020,27246,27466,27686,27898,28106,28312,28512,28708,
  28898,29084,29268,29448,29622,29792,29956,30118,30274,30426,30572,30714,30852,30986,
  31114,31238,31356,31470,31582,31686,31784,31880,31972,32056,32138,32214,32286,32352,
  32414,32470,32522,32568,32610,32646,32680,32706,32728,32746,32758,32764,32767,32764,
  32758,32746,32728,32706,32680,32646,32610,32568,32522,32470,32414,32352,32286,32214,
  32138,32056,31972,31880,31784,31686,31582,31470,31356,31238,31114,30986,30852,30714,
  30572,30426,30274,30118,29956,29792,29622,29448,29268,29084,28898,28708,28512,28312,
  28106,27898,27686,27466,27246,27020,26792,26558,26320,26076,25832,25582,25330,25074,
  24812,24546,24282,24010,23730,23452,23170,22886,22594,22302,22006,21706,21404,21096,
  20788,20476,20158,19842,19520,19196,18868,18536,18206,17868,17530,17190,16846,16498,
  16152,15800,15446,15090,14732,14372,14012,13644,13278,12910,12540,12166,11794,11416,
  11040,10660,10280,9896,9512,9126,8740,8352,7962,7572,7180,6786,6394,5996,5604,5206,
  4808,4410,4010,3612,3212,2812,2412,2008,1608,1206,802,404};


#if LIIROPT
/**
 * @brief 为高通滤波器分配内存
 * 
 * @param HPF_inst &(void* inst)
 * @return int 分配失败返回-1, 否则返回0
 */
int InnoTalkHpf_Create(void **HPF_inst)
{
  nds_liir_f32_t *stt;
  if(HPF_inst == NULL)
  {
    return -1;
  }
  stt = (nds_liir_f32_t *)malloc(sizeof(nds_liir_f32_t));
  *HPF_inst = stt;
  if(stt == NULL)
  {
    return -1;
  }
  return 0;
}

/**
 * @brief 初始化高通滤波器
 * 
 * @param inst1 指向高通滤波器的指针
 * @param rcoeff 高通滤波器的r系数, 和matlab的k系数倒序
 * @param lcoeff 高通滤波器的l系数, 和matlab的v系数倒序
 * @param nstage 高通滤波器阶数
 * @return int32_t 初始化失败则返回-1, 否则返回0
 */
int32_t InnoTalkHpf_InitCore(void* inst1, float* rcoeff, float* lcoeff, int nstage)
{
  nds_liir_f32_t* inst = (nds_liir_f32_t *)inst1;
  if(inst == NULL)
  {
    return -1;
  }
  inst->nstage = nstage;
  inst->state = (float*)malloc(sizeof(float) * (nstage + FRAME_LEN));
  inst->rcoeff = (float*)malloc(sizeof(float) * nstage);
  inst->lcoeff = (float*)malloc(sizeof(float) * (nstage + 1));
  memcpy(inst->rcoeff, rcoeff, nstage * sizeof(float));
  memcpy(inst->lcoeff, lcoeff, (nstage + 1) * sizeof(float));
  memset(inst->state, 0, (nstage + FRAME_LEN) * sizeof(float));

  return 0;

}

/**
 * @brief 对输入帧做高通滤波
 * 
 * @param inst1 指向高通滤波器的指针
 * @param inFrame 输入的一帧short型数据, 帧长128
 * @param outFrame 输出的一帧short型数据, 帧长128
 * @return int 滤波成功返回0
 */
int InnoTalkHpf_ProcessFloat(void* inst1, short* inFrame, short* outFrame)
{
  nds_liir_f32_t* inst = (nds_liir_f32_t *)inst1;
  int i;
  float flInL[FRAME_LEN] = { 0.0 };
  float flOutL[FRAME_LEN] = { 0.0 };

  for (i = 0; i < FRAME_LEN; i++)
	{
		flInL[i] = (float)inFrame[i] / 0x8000;
	}
  nds_liir_f32(inst, flInL, flOutL, FRAME_LEN);
  for (i = 0; i < FRAME_LEN; i++)
  {
    outFrame[i] = (short)(32768 * flOutL[i] + 0.5);
  }
  return 0;
}
#endif

#if FIRAVERAGE
/**
 * @brief 为FIR平均滤波器分配内存
 * 
 * @param smooth_inst &(void* inst)
 * @return int 分配失败返回-1, 否则返回0
 */
int InnoTalkSmooth_Create(void **smooth_inst)
{
  nds_fir_q15_t *stt;
  if(smooth_inst == NULL)
  {
    return -1;
  }
  stt = (nds_fir_q15_t *)malloc(sizeof(nds_fir_q15_t));
  *smooth_inst = stt;
  if(stt == NULL)
  {
    return -1;
  }
  return 0;
}

/**
 * @brief 初始化一个FIR平均滤波器
 * 
 * @param inst1 FIR滤波器指针
 * @param coeff_size 系数的长度, 即窗长
 * @return int32_t 初始化失败返回-1, 否则返回0
 */
int32_t InnoTalkSmooth_InitCore(void* inst1, int coeff_size)
{
  nds_fir_q15_t* inst = (nds_fir_q15_t *)inst1;
  if(inst == NULL)
  {
    return -1;
  }
  int16_t val;
  val = (int16_t)((int32_t)(Q15MOD) / (int32_t)(coeff_size));
  inst->coeff_size = (uint32_t)coeff_size;
  inst->state = (int16_t*)malloc(sizeof(int16_t) * (coeff_size + HALF_ANAL_BLOCKL - 1));
  inst->coeff = (int16_t*)malloc(sizeof(int16_t) * coeff_size);
  memset(inst->state, 0, (coeff_size + HALF_ANAL_BLOCKL - 1) * sizeof(int16_t));
  memset(inst->coeff, val, coeff_size * sizeof(int16_t));

  return 0;

}
#endif

#if FSmooth
/**
 * @brief 创建一个滑动平均滤波器并初始化
 * 
 * @param size 滑动窗口长度
 * @return MovingAverage* 滑动平均滤波器obj
 */
MovingAverage* movingAverageCreate(int32_t size) {
    MovingAverage* obj = (MovingAverage*)malloc(sizeof(MovingAverage));
    obj->size = size;
    obj->sum = 0;
    obj->queue = (int16_t*)malloc(sizeof(int16_t) * (size + 1));
    obj->front = 0;
    obj->rear = 0;
    return obj;
}

/**
 * @brief 对输入的值做滑动平均
 * 
 * @param obj 滑动平均滤波器
 * @param val 输入的下一个值
 * @return int16_t 滑动平均后的值
 */
int16_t movingAverageNext(MovingAverage* obj, int16_t val) {
    int size = (obj->rear - obj->front + obj->size + 1) % (obj->size + 1);
    if (size == obj->size) {
        obj->sum -= obj->queue[obj->front];
        obj->front = (obj->front + 1) % (obj->size + 1);
        size--;
    }
    obj->queue[obj->rear] = val;
    obj->rear = (obj->rear + 1) % (obj->size + 1);
    obj->sum += val;
    size++;
    return (int16_t)(obj->sum / size);
}

/**
 * @brief 释放滑动平均滤波器的内存空间
 * 
 * @param obj 滑动平均滤波器
 */
void movingAverageFree(MovingAverage* obj) {
    free(obj->queue);
    free(obj);
}
#endif

/**
 * @brief 为ns结构体分配内存空间
 * 
 * @param NS_inst &(void *)
 * @return int 分配失败返回-1, 否则返回0
 */
int InnoTalkNsx_Create(void **NS_inst)
{
	NsxInst_t *stt;
	if (NS_inst == NULL)
	{
		return -1;
	}
	stt = (NsxInst_t *)malloc(sizeof(NsxInst_t));

	*NS_inst = stt;
	if (stt == NULL)
	{
		return -1;
	}
	return 0;
}

/**
 * @brief 初始化ns结构体
 * 
 * @param inst1 指向ns结构体的指针
 * @param fs 采样频率
 * @return int32_t 初始化失败返回-1, 否则返回0
 */
int32_t InnoTalkNsx_InitCore(void* inst1, uint32_t fs)
{
  NsxInst_t* inst = (NsxInst_t *)inst1;
  if (inst == NULL) {
    return -1;
  }
  // Initialization of struct
  inst->fs = fs;
  inst->blockLen = 128;
  inst->stages = 8;
  inst->blockIndex = -1; //frame counter
  inst->scaleEnergyIn = 0;
  inst->energyIn = 0;
  inst->zeroInputSignal = 0;
  inst->norm = 0;
  inst->normPrev = 15;
  inst->overdrive16 = 42598;   // Q15(1.3)
  inst->denoiseBound16 = 2294; // Q15(0.07)

  nds_set_q15((int16_t)0, inst->analysisBuffer, (uint32_t)ANAL_BLOCKL_MAX);
  nds_set_q15((int16_t)0, inst->synthesisBuffer, (uint32_t)ANAL_BLOCKL_MAX);
  nds_set_q15((int16_t)0, inst->fftdata, (uint32_t)(ANAL_BLOCKL_MAX * 2));
  nds_set_q15((int16_t)0, inst->smooth32, (uint32_t)HALF_ANAL_BLOCKL);
  nds_set_q31((int32_t)0, inst->pmagnPrev32, (uint32_t)HALF_ANAL_BLOCKL);
  nds_set_q31((int32_t)0, inst->noisePrev32, (uint32_t)HALF_ANAL_BLOCKL);
  nds_set_q31((int32_t)0, inst->probPrev32, (uint32_t)HALF_ANAL_BLOCKL);
  nds_set_q31((int32_t)0, inst->signalPrev32, (uint32_t)HALF_ANAL_BLOCKL);
  nds_set_q31((int32_t)0, inst->minMagn32, (uint32_t)HALF_ANAL_BLOCKL);

  inst->initFlag = 1;
  return 0;
}

#if CUT
int InnoTalkSpl_GetScalingSquare(int16_t *in_vector, int in_vector_length, int times)
{
	int nbits = InnoTalkSpl_GetSizeInBits(times);
	int i;
	int16_t smax = -1;
	int16_t sabs;
	int16_t *sptr = in_vector;
	int t;
	int looptimes = in_vector_length;

	for (i = looptimes; i > 0; i--)
	{
		sabs = (*sptr > 0 ? *sptr++ : -*sptr++);
		smax = (sabs > smax ? sabs : smax);
	}
	t = InnoTalkSpl_NormW32(INNOTALK_SPL_MUL(smax, smax));

	if (smax == 0)
	{
		return 0; // Since norm(0) returns 0
	}
	else
	{
		return (t > nbits) ? 0 : nbits - t;
	}
}
#endif

/**
 * @brief 计算信号能量
 * 
 * @param vector 输入int16的信号
 * @param vector_length 信号长度
 * @param scale_factor 信号归一化尺度
 * @return uint64_t 信号能量
 */
uint64_t InnoTalk_Energy(int16_t* vector, uint32_t vector_length, int* scale_factor)
{
	uint64_t en = 0;

#if CUT
	int scaling =  InnoTalkSpl_GetScalingSquare(vector, vector_length, vector_length);
#else
	int scaling = 0;
#endif

	en = (uint64_t)nds_pwr_q15(vector,vector_length);
	en = en >> scaling;
	*scale_factor = scaling;
	return en;
}

/**
 * @brief 对输入信号进行加窗, 缓存, 计算能量, 移位以及FFT
 * 
 * @param inst1 ns结构体的指针
 * @param speechFrame 输入的一帧信号
 */
void InnoTalkNsx_DataAnalysis(void *inst1, short* speechFrame)
{
	NsxInst_t* inst = (NsxInst_t *)inst1;
	int i;
  int16_t  maxWinData;
  int16_t  winData[ANAL_BLOCKL_MAX]= { 0 }, realImag[ANAL_BLOCKL_MAX] = { 0 };
  int16_t  DataAbs[ANAL_BLOCKL_MAX] = { 0 };
  uint32_t dataindex = 0;

  nds_dup_q15(inst->analysisBuffer + inst->blockLen, inst->analysisBuffer, (uint32_t)(ANAL_BLOCKL_MAX - inst->blockLen));
  nds_dup_q15(speechFrame, inst->analysisBuffer + ANAL_BLOCKL_MAX - inst->blockLen, (uint32_t)inst->blockLen);

  nds_mul_q15(kBlocks128w256x, inst->analysisBuffer, winData, (uint32_t)ANAL_BLOCKL_MAX);
  // nds_shift_q15(winData,1,winData,ANAL_BLOCKL_MAX);
#if EngScaled
  // Get input energy
  inst->energyIn = InnoTalk_Energy(winData, (uint32_t)ANAL_BLOCKL_MAX, &(inst->scaleEnergyIn));
#endif
  // Reset zero input flag
  inst->zeroInputSignal = 0;
  // Acquire norm for winData
  nds_abs_q15(winData, DataAbs, (uint32_t)ANAL_BLOCKL_MAX);
  maxWinData = nds_max_q15(DataAbs, (uint32_t)ANAL_BLOCKL_MAX, &dataindex);
  inst->norm = InnoTalkSpl_NormW16(maxWinData);
  if (maxWinData == 0) {
    // Treat zero input separately.
    inst->zeroInputSignal = 1;
    return;
  }

  nds_shift_q15(winData, (int8_t)inst->norm, realImag, (uint32_t)ANAL_BLOCKL_MAX);

  for (i = 0; i < ANAL_BLOCKL_MAX; i++)
  {
	  FFTIOtT[2 * i] = realImag[i];
	  FFTIOtT[2 * i + 1] = 0;
  }
  nds_cfft_rd4_q15(FFTIOtT, inst->stages);
  nds_dup_q15(FFTIOtT, inst->fftdata, (uint32_t)(2 * ANAL_BLOCKL_MAX));
  // nds_shift_q15(FFTIOtT, 1, inst->fftdata, 2 * ANAL_BLOCKL_MAX);
}

/**
 * @brief 降噪算法的核心函数
 * 
 * @param inst1 ns结构体
 * @param speechFrame 输入的short型一帧信号
 * @param outFrame 输出的short型一帧降噪后信号
 * @return int 降噪失败返回-1, 否则返回0
 */
int InnoTalkNsx_ProcessCore(void* inst1, short* speechFrame, short* outFrame)
{
  // main routine for noise suppression
  NsxInst_t* inst = (NsxInst_t *)inst1;
  if (inst->initFlag != 1) {
    return -1;
  }

  int i;
  int16_t  dTmp, Qdiff;
  int32_t  fTmp32;
  int64_t  tempdata64;
  int32_t  vecTmp32_1[HALF_ANAL_BLOCKL], vecTmp32_2[HALF_ANAL_BLOCKL];
  int16_t  vecTmp16_1[HALF_ANAL_BLOCKL], vecTmp16_2[HALF_ANAL_BLOCKL];

  // int32_t  winDataI[HALF_ANAL_BLOCKL * 2];
  // int32_t  winDataTmp[ANAL_BLOCKL_MAX];
  int16_t  winDataI[ANAL_BLOCKL_MAX * 2];
  int16_t  winDataO[ANAL_BLOCKL_MAX];
  int16_t  real16[HALF_ANAL_BLOCKL],imag16[HALF_ANAL_BLOCKL];
  int32_t  magn32[HALF_ANAL_BLOCKL], noise32[HALF_ANAL_BLOCKL], pmagn32[HALF_ANAL_BLOCKL];
  int64_t  snrLocPost32[HALF_ANAL_BLOCKL] = { 0 };
  int64_t  snrLocPrior32[HALF_ANAL_BLOCKL] = { 0 };
  int32_t  probSpeechFinal32[HALF_ANAL_BLOCKL] = { 0 }; /* Q15 */
  int32_t  ii32[HALF_ANAL_BLOCKL] = {0}; /* Q15 */
  int16_t  theFilter32[HALF_ANAL_BLOCKL];
  int32_t  pr32[HALF_ANAL_BLOCKL] = {0}, afa32[HALF_ANAL_BLOCKL] = { 0 }; /* Q15 */
  int32_t  factor32, gain32;
  uint64_t energyOut;
#if FSmooth
  uint64_t fenergy1, fenergy2;
  int32_t  zeta;
  int32_t  winLen;
  MovingAverage* smooth;
#endif
#if (FIRAVERAGE)
  void * smooth_inst = NULL;
#endif

  // Store speechFrame and transform to frequency domain
  InnoTalkNsx_DataAnalysis(inst, speechFrame);
  Qdiff = 2 * (inst->norm - inst->normPrev);

  if (inst->zeroInputSignal) {
    // synthesize the special case of zero input
    // read out fully processed segment
    nds_dup_q15(inst->synthesisBuffer, outFrame, (uint32_t)inst->blockLen);
    // update synthesis buffer
    nds_dup_q15(inst->synthesisBuffer + inst->blockLen, inst->synthesisBuffer, (uint32_t)(ANAL_BLOCKL_MAX - inst->blockLen));
    nds_set_q15((int16_t)0, inst->synthesisBuffer + ANAL_BLOCKL_MAX - inst->blockLen, (uint32_t)inst->blockLen);
    return 0;
  }

  inst->blockIndex++; // Update the block index only when we process a block.

  // Q(norm)
  imag16[0] = inst->fftdata[1];
  real16[0] = inst->fftdata[0];
  imag16[HALF_ANAL_BLOCKL - 1] = inst->fftdata[ANAL_BLOCKL_MAX + 1];
  real16[HALF_ANAL_BLOCKL - 1] = inst->fftdata[ANAL_BLOCKL_MAX];
  // Q(2 * norm)
  magn32[0] = (int32_t)(real16[0] * real16[0]);
  magn32[HALF_ANAL_BLOCKL - 1] = (int32_t)(real16[HALF_ANAL_BLOCKL - 1] * real16[HALF_ANAL_BLOCKL - 1]);
  // 计算输入的频域能量, Q(2 * norm)
#if FSmooth
  fenergy1 = (uint64_t)magn32[0];
  fenergy1 += (uint64_t)magn32[HALF_ANAL_BLOCKL - 1];
#endif

  for (i = 1; i < HALF_ANAL_BLOCKL - 1; i++)
  {
    // Q(norm)
    real16[i] = inst->fftdata[2 * i];
    imag16[i] = inst->fftdata[2 * i + 1];
    // Q(2 * norm)
    fTmp32 = (int32_t)(real16[i] * real16[i]);
    fTmp32 += (int32_t)(imag16[i] * imag16[i]);
    magn32[i] = fTmp32;
#if FSmooth
    fenergy1 += (uint64_t)magn32[i];
#endif
  }

  // 如果是第一帧
  if (inst->blockIndex == 0)
  {
    nds_dup_q31(magn32, pmagn32, (uint32_t)HALF_ANAL_BLOCKL);
    nds_dup_q31(magn32, inst->minMagn32, (uint32_t)HALF_ANAL_BLOCKL);
    nds_dup_q31(magn32, noise32, (uint32_t)HALF_ANAL_BLOCKL);
    nds_set_q31((int32_t)0, probSpeechFinal32, (uint32_t)HALF_ANAL_BLOCKL);
  }
  else
  // 非首帧
  {
    // norm < normPrev, 当前帧的magn右移位数更多, 则把上一帧数据再右移
    if (Qdiff < 0)
    {
      nds_shift_q31(inst->pmagnPrev32, (int8_t)Qdiff, inst->pmagnPrev32, (uint32_t)HALF_ANAL_BLOCKL); //Q(2*norm)
      nds_shift_q31(inst->minMagn32, (int8_t)Qdiff, inst->minMagn32, (uint32_t)HALF_ANAL_BLOCKL); //Q(2*norm)
      nds_shift_q31(inst->noisePrev32, (int8_t)Qdiff, inst->noisePrev32, (uint32_t)HALF_ANAL_BLOCKL); //Q(2*norm)
      nds_shift_q31(inst->signalPrev32, (int8_t)Qdiff, inst->signalPrev32, (uint32_t)HALF_ANAL_BLOCKL); //Q(2*norm)
    }
    // norm > normPrev, 上一帧帧的magn右移位数更多, 则把当前帧数据再右移
    else if (Qdiff > 0) // norm > normPrev, 转到normPrev
    {
      nds_shift_q31(magn32, (int8_t)(-Qdiff), magn32, (uint32_t)HALF_ANAL_BLOCKL);
    }

    // 1. 平滑功率谱, 对应原理中最小值统计估噪的第(1)步即公式(1)
    nds_scale_q31(inst->pmagnPrev32, (int32_t)SMOOTH_APY16, (int8_t)16, vecTmp32_1, (uint32_t)HALF_ANAL_BLOCKL);
    nds_scale_q31(magn32, (int32_t)SMOOTH_APY16_S, (int8_t)16, vecTmp32_2, (uint32_t)HALF_ANAL_BLOCKL);
    nds_add_q31(vecTmp32_1, vecTmp32_2, pmagn32, (uint32_t)HALF_ANAL_BLOCKL);
 	  for (i = 0; i < HALF_ANAL_BLOCKL; i++)
  	{
  	  // // 1.平滑功率谱 对应原理中最小值统计估噪的第(1)步即公式(1)
  	  // tempdata64 = (int64_t)(SMOOTH_APY16 * (int64_t)inst->pmagnPrev32[i]) + (int64_t)(SMOOTH_APY16_S * (int64_t)magn32[i]); //Q(min + 15)
  	  // pmagn32[i] = (int32_t)(tempdata64 >> 15); //Q(min)


      // 2.搜索频带最小值 对应原理中最小值统计估噪的第(2)步
  	  if (inst->minMagn32[i] < pmagn32[i])
  	  {
        tempdata64 = (int64_t)(SMOOTH_R16 * (int64_t)inst->minMagn32[i]) + (int64_t)(1638 * ((int64_t)pmagn32[i] - (((int64_t)(SMOOTH_BETA16 * (int64_t)inst->pmagnPrev32[i])) >> 15))); //Q(min + 15)
        inst->minMagn32[i] = (int32_t)(tempdata64 >> 15); //Q(min)
      }
  	  else
  		{
  	    inst->minMagn32[i] = pmagn32[i]; //Q(min)
  		}

      // 3.判断是否存在语音 对应原理中最小值统计估噪的第(3)步
  		if (inst->minMagn32[i] <= 0)
  		{
  		  pr32[i] = delta32[i] + Q15MOD;
  		}
  		else
  		{
  		  pr32[i] = (((int64_t)pmagn32[i]) << 15) / inst->minMagn32[i];
  		}
  		if (pr32[i] > delta32[i])
  		{
  		  ii32[i] = (int32_t)Q15MOD;
  		}
  		else
  		{
  		  ii32[i] = 0;
  		}

  		// // 4.计算语音出现的概率 对应原理中最小值统计估噪的第(4)步
      // tempdata64 = (int64_t)AP16 * inst->probPrev32[i] + (int32_t)(AP16_S) * ii32[i]; // Q15 * Q15
      // probSpeechFinal32[i] = (int32_t)(tempdata64 >> 15); // Q15
  	}

    // 4.计算语音出现的概率 对应原理中最小值统计估噪的第(4)步
    nds_scale_q31(inst->probPrev32, (int32_t)AP16, (int8_t)16, vecTmp32_1, (uint32_t)HALF_ANAL_BLOCKL);
    nds_scale_q31(ii32, (int32_t)AP16_S, (int8_t)16, vecTmp32_2, (uint32_t)HALF_ANAL_BLOCKL);
    nds_add_q31(vecTmp32_1, vecTmp32_2, probSpeechFinal32, (uint32_t)HALF_ANAL_BLOCKL);

    // 5.噪声谱估计,对应原理中最小值统计估噪的第(5)步:
    if (inst->blockIndex < 7)
    {
      nds_dup_q31(magn32, noise32, (uint32_t)HALF_ANAL_BLOCKL);
    }
    else
    {
      for (i = 0; i < HALF_ANAL_BLOCKL; i++)
      {
        tempdata64 = (int64_t)G16 * probSpeechFinal32[i] + (int64_t)966367641;//(0.9) * Q30MOD
        afa32[i] = (int32_t)(tempdata64 >> 15); //Q15
        tempdata64 = (int64_t)(afa32[i] * (int64_t)inst->noisePrev32[i]) + (int64_t)(Q15MOD - afa32[i]) * (int64_t)magn32[i]; // Q(min + 15)
        noise32[i] = (int32_t)(tempdata64 >> 15);  //Q(min)
      }
      // nds_scale_q31(probSpeechFinal32, (int32_t)G16, (int8_t)16, Vectmp1, (uint32_t)HALF_ANAL_BLOCKL);
      // nds_offset_q31(Vectmp1, (int32_t)G16_S, afa32, (uint32_t)HALF_ANAL_BLOCKL);
      // nds_neg_q31(afa32, Vectmp2, (uint32_t)HALF_ANAL_BLOCKL);
      // nds_offset_q31(Vectmp2, (int32_t)Q15MOD, Vectmp2, (uint32_t)HALF_ANAL_BLOCKL);
      // nds_mul_q31(afa32, inst->noisePrev32, Vectmp1, (uint32_t)HALF_ANAL_BLOCKL);
      // nds_mul_q31(Vectmp2, magn32, Vectmp2, (uint32_t)HALF_ANAL_BLOCKL);
      // nds_add_q31(Vectmp1, Vectmp2, noise32, (uint32_t)HALF_ANAL_BLOCKL);
    }
  }

  // 将当前帧的值赋给前一帧
  nds_dup_q31(probSpeechFinal32, inst->probPrev32, (uint32_t)HALF_ANAL_BLOCKL);
  if (Qdiff < 0)
  {
    nds_dup_q31(pmagn32, inst->pmagnPrev32, (uint32_t)HALF_ANAL_BLOCKL);
    nds_dup_q31(noise32, inst->noisePrev32, (uint32_t)HALF_ANAL_BLOCKL);
  }
  else
  {
    nds_shift_q31(pmagn32, (int8_t)Qdiff, inst->pmagnPrev32, (uint32_t)HALF_ANAL_BLOCKL);
    nds_shift_q31(noise32, (int8_t)Qdiff, inst->noisePrev32,(uint32_t)HALF_ANAL_BLOCKL);
    nds_shift_q31(inst->minMagn32, (int8_t)Qdiff, inst->minMagn32, (uint32_t)HALF_ANAL_BLOCKL);
  }
  inst->normPrev = inst->norm;

  // 6. 计算先验/后验信噪比和维纳滤波器
  if (inst->blockIndex == 0)
  {
    nds_set_q15((int16_t)13919, theFilter32, (uint32_t)HALF_ANAL_BLOCKL);
    nds_scale_q31(magn32, (int32_t)5912,(int8_t)16, inst->signalPrev32, (uint32_t)HALF_ANAL_BLOCKL);
  }
  else
  {
    for ( i = 0; i < HALF_ANAL_BLOCKL; i++)
    {
      if (noise32[i] == 0)
      {
        theFilter32[i] = (int16_t)13757;
      }
      else
      {
        snrLocPost32[i] = magn32[i]; //Q(min)
        tempdata64 = (int64_t)(AF16_S) * INNOTALK_SPL_MAX((int32_t)snrLocPost32[i] - (int32_t)noise32[i], 0); //Q(min + 15)
        snrLocPrior32[i] = (int64_t)(AF16 * (int64_t)inst->signalPrev32[i]) + tempdata64; //Q(min + 15)
        tempdata64 = (int64_t)inst->overdrive16 * noise32[i]; //Q(min + 15)
        theFilter32[i] = ((int64_t)snrLocPrior32[i] * Q15MOD) / (tempdata64 + snrLocPrior32[i]); //Q((min + 30) - (min + 15)) = Q15
        if (theFilter32[i] > TMPVALUE)
        {
          theFilter32[i] = TMPVALUE;
        }
        if (theFilter32[i] < inst->denoiseBound16)
        {
          theFilter32[i] = inst->denoiseBound16;
        }
      }
      tempdata64 = ((int64_t)(((int64_t)magn32[i] * (int64_t)theFilter32[i]) >> 15)); // Q(min + 15 - 15)
      inst->signalPrev32[i] = (int32_t)(((int64_t)theFilter32[i] * tempdata64) >> 15);  // Q(min + 15 - 15)
    }
  }

  if (Qdiff > 0)
  {
  	nds_shift_q31(inst->signalPrev32, (int8_t)Qdiff, inst->signalPrev32, (uint32_t)HALF_ANAL_BLOCKL);
  }

  // 7. 抑制水声
  nds_scale_q15(inst->smooth32, (int16_t)FILTER_SMOOTH, (int8_t)0, vecTmp16_1, (uint32_t)HALF_ANAL_BLOCKL);
  nds_scale_q15(theFilter32, (int16_t)FILTER_SMOOTH_S, (int8_t)0, vecTmp16_2, (uint32_t)HALF_ANAL_BLOCKL);
  nds_add_q15(vecTmp16_1, vecTmp16_2, inst->smooth32, (uint32_t)HALF_ANAL_BLOCKL);
#if FSmooth
  fenergy2 = 0;
  for (i = 0; i < HALF_ANAL_BLOCKL; i++)
  {
    fenergy2 += (uint64_t)inst->signalPrev32[i];
  }
  zeta = (fenergy2 << 15) / fenergy1;
  if (zeta > ZETA_THR)
  {
    zeta = (int32_t)Q15MOD;
  }
  if (zeta == (int32_t)Q15MOD)
  {
    winLen = 1;
  }
  else
  {
    winLen = 2 * ((((Q15MOD - (zeta << 15) / ZETA_THR) * PSI) + Q15MOD) >> 15) + 1;
    smooth = movingAverageCreate(winLen);
    for (i = 0; i < HALF_ANAL_BLOCKL; i++)
    {
      theFilter32[i] = movingAverageNext(smooth, theFilter32[i]);
    }
    movingAverageFree(smooth);
    // InnoTalkSmooth_Create(&smooth_inst);
    // InnoTalkSmooth_InitCore(smooth_inst, winLen);
    // nds_fir_q15(smooth_inst, theFilter32, theFilter32, HALF_ANAL_BLOCKL);
  }
  memcpy(inst->smooth32, theFilter32, HALF_ANAL_BLOCKL * sizeof(int16_t));
#endif

  // 8. 增强语音谱
  nds_mul_q15(inst->smooth32, real16, real16, (uint32_t)HALF_ANAL_BLOCKL);
  nds_mul_q15(inst->smooth32, imag16, imag16, (uint32_t)HALF_ANAL_BLOCKL);

  winDataI[0] = real16[0];
  winDataI[1] = imag16[0];
  winDataI[ANAL_BLOCKL_MAX] = real16[HALF_ANAL_BLOCKL - 1];
  winDataI[ANAL_BLOCKL_MAX + 1] = imag16[HALF_ANAL_BLOCKL - 1];
  for (i = 1; i < HALF_ANAL_BLOCKL - 1; i++) {
	  winDataI[2 * i] = real16[i];
	  winDataI[2 * i + 1] = imag16[i];
  }

  for(i = HALF_ANAL_BLOCKL - 2;i>0;i--)
  {
	  winDataI[ANAL_BLOCKL_MAX * 2 - i * 2] = real16[i];
	  winDataI[ANAL_BLOCKL_MAX * 2 - i * 2 + 1] = -imag16[i];
  }

  nds_cifft_rd4_q15(winDataI,inst->stages);
  for (i = 0; i < ANAL_BLOCKL_MAX; i++)
  {
	  winDataO[i] = winDataI[i * 2];
  }
  nds_shift_q15(winDataO, (int8_t)(-inst->norm), winDataO, (uint32_t)ANAL_BLOCKL_MAX);

  factor32 = Q15MOD;
  if (inst->blockIndex > 0) {
	  energyOut = (uint64_t)nds_pwr_q15(winDataO, (uint32_t)ANAL_BLOCKL_MAX);
	  energyOut = (uint64_t)(energyOut << 30);
	  fTmp32 = (int32_t)(energyOut / (inst->energyIn + 32768));
	  gain32 = nds_sqrt_f32((float)fTmp32);
	  if (gain32 > B_LIM16) {
	    factor32 = (int32_t)Q15MOD + ((42598 * (gain32 - B_LIM16)) >> 15);
	  }
	  else {
	    factor32 = (int32_t)Q15MOD + 2 * (gain32 - B_LIM16);
	  }
	  if (gain32 * factor32 > (int64_t)Q30MOD) {
	    factor32 = (int64_t)Q30MOD / gain32;
	  }
  }


  for (i = 0; i < ANAL_BLOCKL_MAX; i++) {
	  inst->synthesisBuffer[i] += (short)(((int64_t)factor32*(int32_t)kBlocks128w256x[i] * (int32_t)winDataO[i]) >> 30);
  }

  for (i = 0; i < inst->blockLen; i++) {
	  dTmp = (short)inst->synthesisBuffer[i];
	  outFrame[i] = dTmp;
  }
  // update synthesis buffer
  nds_dup_q15(inst->synthesisBuffer + inst->blockLen, inst->synthesisBuffer, (uint32_t)(ANAL_BLOCKL_MAX - inst->blockLen));
  nds_set_q15((int16_t)0, inst->synthesisBuffer + ANAL_BLOCKL_MAX - inst->blockLen, (uint32_t)inst->blockLen);

  return 0;
}
