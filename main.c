/**
 * @file main.c
 * @author fjj
 * @brief Telink对语音信号做处理的仿真程序
 * @version 0.1
 * @date 2022-08-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "nsx.h"
// #include "agc.h"


#if (LIIROPT)
#define nstage        (3)   //hpf阶数
#endif

/**
 * @brief 对输入音频做算法处理的函数
 * 
 * @param szFileIn 输入的文件名
 * @param szFileOut 输出的文件名
 * @param nSample 采样率
 */
void InnoTalkNS16KSampleX(char *szFileIn, char *szFileOut, uint32_t nSample)
{

#if (LIIROPT)
	void *HPF_inst = NULL;
#endif
	void *pNS_inst = NULL;
	// void *AGC_inst = NULL;
	short shInL[FRAME_LEN] = {0};
	// short shTmpL[FRAME_LEN] = { 0 };
	short shOutL[FRAME_LEN] = {0};
	FILE *fpIn = NULL;
	FILE *fpOut = NULL;

	do
	{
		fpIn = fopen(szFileIn, "rb");
		fpOut = fopen(szFileOut, "wb");
		if (NULL == fpIn || NULL == fpOut)
		{
			printf("open file err \n");
			break;
		}
#if (LIIROPT)
		float rcoeff[nstage] = {-0.8876, 0.9983, -0.9989}; // hpf r系数, 和matlab的k系数倒序
		float lcoeff[nstage + 1] = {-0.94252, 0.11082, 0.00464, -0.000016}; //hpf l系数, 和matlab的v系数倒序
		InnoTalkHpf_Create(&HPF_inst);
		InnoTalkHpf_InitCore(HPF_inst, rcoeff, lcoeff, nstage);
#endif
		InnoTalkNsx_Create(&pNS_inst);
		InnoTalkNsx_InitCore(pNS_inst, nSample);

		// InnoTalkAgc_Create(&AGC_inst);
		// int32_t minLevel = 0;
		// int32_t maxLevel = 255;
		// // int agcMode = kAgcModeAdaptiveDigital;
		// InnoTalkAgc_Init(AGC_inst, minLevel, maxLevel, nSample);
		// InnoTalkAgc_config_t agcConfig;
		// agcConfig.compressionGaindB = (int16_t)15;
		// agcConfig.limiterEnable = (uint8_t)1;
		// agcConfig.targetLevelDbfs = (int16_t)6;
		// agcConfig.SilenceGainFall = (int16_t)110;
		// agcConfig.AgcVadLowThr = (int16_t)0;                //JT:AGC中VAD判定低门限,Q10，推荐值为0
		// agcConfig.AgcVadUppThr = (int16_t)1024;             //JT:AGC中VAD判定高门限,Q10，推荐值为1024
		// agcConfig.GateLowThr = (int16_t)0;                  //JT:Gate判定低门限，推荐值为0
		// agcConfig.GateUppThr = (int16_t)2500;               //JT:Gate判定高门限，推荐值为2500
		// agcConfig.GateVadSensitivity = (int16_t)6;             //JT:Gate计算对语音活性的敏感度，推荐值为6
		// InnoTalkAgc_set_config(AGC_inst, agcConfig);

		while (1)
		{
			if (FRAME_LEN == fread(shInL, sizeof(short), FRAME_LEN, fpIn))
			{
				InnoTalkNsx_ProcessCore(pNS_inst, shInL, shOutL);
				// InnoTalkAgc_Process(AGC_inst, shTmpL, shOutL, agcConfig);
				fwrite(shOutL, sizeof(short), FRAME_LEN, fpOut);
			}
			else
			{
				break;
			}
		}
	} while (0);

	fclose(fpIn);
	fclose(fpOut);
}


int main()
{
	printf("processing...\n");
	InnoTalkNS16KSampleX("TelinkTest.pcm", "TelinkTest_0920.pcm", (uint32_t)16000);
	printf("end!\n");

	getchar();
	return 0;
}
