#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES	// pi값 사용. 변수 이름 M_PI
#include <math.h>
#include <errno.h>
#include <direct.h>
//엔탈피나 cp 값 kj에서 j 단위로 다 치환할 것.
#define GRID 100 //	쪼개놓은 Partition을 기준으로 온도 테이블을 GRID 수만큼 더 쪼개서 테이블로 기입.
//GRID 값이 100보다 작은 범위에서 Hot_in_temp과 Cold_out_temp가 겹쳐서 에러를 일으킬 가능성이 있다.
//CSV 파일로 저장할 때 geometric이나 그냥 thermal condition들은 일단 txt파일로 다 옮겨놓을 것.
//대부분 데이터들을 double로 처리할 것
//<Database>
	//<Struction>
typedef struct _tempdata{ // ProData
	double temperature;
	double conductivity;
	double viscosity;
	double density;
	double Cp;
	double enthalpy;
}ProData;	//Property_Data

typedef struct _hotside{ // Tube
	//온도 관련 프로퍼티
	double in_temp;
	double out_temp;
	double avg_temp;
	int index;	// 온도 테이블에 대한 주소값.
	double in_enthalpy;
	double out_enthalpy;// 엔탈피만 인 아웃을 구분하고 나머지는 다 평균 온도 기반으로 계산.
	double conductivity;
	double viscosity;
	double density;
	double Cp;
	//그 외 프로퍼티
	double Pr;
	double e_ratio;
	double p_ratio;
	double angle_ratio;
	double velocity;
	double Re;
	double Nu;
	double hi;
	double fi;
	double DP;
}Tube;

typedef struct _coldside{ //Annular
	//온도 관련 프로퍼티
	double in_temp;
	double out_temp;
	double avg_temp;
	int index;	// 온도 테이블에 대한 주소값.
	double in_enthalpy;
	double out_enthalpy;// 엔탈피만 인 아웃을 구분하고 나머지는 다 평균 온도 기반으로 계산.
	double conductivity;
	double viscosity;
	double density;
	double Cp;
	//그 외 프로퍼티
	double Pr;
	double e_ratio;
	double p_ratio;
	double angle_ratio;
	double r_ratio;
	double Dhyd;
	double Aeff;
	double velocity;
	double Re;
	double Nu;
	double ho;
	double ef;
	double fo;
	double DP;
}Annular;

typedef struct _common{ //common
	double HX;
	double length;
	double height;
	double volume;
	double DTln;	//LMTD
	double UA;
	double Ui;
	double Ai;
	double Uo;
	double Ao;
	double Cucond;
}Common;
	//<\Struction>
	//<Parameters>
double T_cond[10] = { 0.5475, 0.00205, -5.55E-06, 4.55E-08, -5.90E-09, 1.56E-10, -2.26E-12, 2.00E-14, -1.00E-16, 2.18E-19 }; // 1st
double Viscosity[10] = { 0.00179, -6.23E-05, 1.70E-06, -3.96E-08, 7.63E-10, -1.14E-11, 1.22E-13, -8.75E-16, 3.69E-18, -6.91E-21 }; // 2nd
double Density[10] = { 999.8292, 0.10526, -0.01371, 3.46E-04, -8.30E-06, 1.28E-07, -1.12E-09, 4.44E-12, 1.17E-15, -4.43E-17 }; // 3rd
double Cp[10] = { 4.22837, -0.0083, 6.28E-04, -2.61E-05, 6.74E-07, -1.15E-08, 1.33E-10, -9.98E-13, 4.35E-15, -8.35E-18 }; // 4th
double Enthalpy[10] = { 0.06374, 4.22664, -0.00381, 1.75E-04, -4.63E-06, 7.13E-08, -6.16E-10, 2.49E-12, -3.65E-16, -2.04E-17 }; // 5th
ProData copper[7];	// 온도-열전달계수만 표기. 나머지 인스턴스들은 건들지 말 것.
	//<\parameters>
	//<Geometric conditions>
//얘네들 값도 따로 저장하는 테이블을 만들어서 저장해서 불러오는 방식도 생각해보자.
double Dbi;	//g
double Deo;	//g
double tw = 5.0E-04;  //constant 
int Ns;	//g
double Pitch;	//g
double Dcan = 0.2;	// constant
double R = 0.6;	//constant
double FB;
double Dvi;
double Dvo;
double Doi;
double angle;
double e;
double mh = 2;
double mc = 1.5;
	//<\Geometric conditions>
	//<Thermal conditions>
double Hot_in_temp;
double Hot_out_temp;
double Cold_in_temp;
double Cold_out_temp;
int partition;	// n, 몇 번을 나눌 것인가??
double DTh;	// 파티션 나눈 갯수에 따라 튜브의 온도 일정한 간격을 담은 값
//double Hot_avg_temp = 57.5;
//double Cold_avg_temp;
	//<\Thermal conditions>
	//<Parametric Analysis Conditions>
double min_factor, max_factor, unit = 0.01;
int min_Ns, max_Ns;
double Dbi_ini, Deo_ini, Doi_ini, Pitch_ini;	// 추후에 위 변수들은 변경될 값들이라 초기값들을 담아주는 변수들.
int Ns_ini; 
int size; // parametric analysis를 위한 크기 저장값.
	//<\Parametric Analysis Conditions>
//<\Database>
/////////////////////////////////////////////////////////////////////////////////////
//<Functions>
	//<Calculator/ShowData/FPrintf Functions>
void write_given_geodata(int option);	//주어지는 값들을 받아서 전역변수에 저장시킨다. 해당 값들은 Geometric 값들, termperature 값들. 따로 체크하는 곳이 필요하다.
void Print_given_geodata(void); // 형상 데이터 확인용 출력
void write_given_tempdata(int option);
void Print_given_tempdata(void);
void write_given_parametric_factor(int option); 	// Parametric Analysis에 쓰일 변수 대입하는 함수.
void Print_given_parametric_factor(int size);
void Initialize_Geoconditions(void);	// Parametric Analysis에 쓰일 Geoconditions을 초기화.
void create_property_csv(ProData* prodata, int add); // 아래 프로퍼티 계산기를 통해 구조체 배열에 저장과 동시에, 프로퍼티 csv파일을 작성한다. 
void property_calculator(ProData* prodata, int add, double temp, double gap);	// 온도에 따른 프로퍼티 값들을 계산한다.
void Calculate_HX_partly(ProData* prodata, Tube* ptube, Annular* pannular, Common* pcommon); //★★★이번 계산기의 총 집합체★★★
void create_counterflow_csv(char* filename, Tube* ptube, Annular* pannular, Common* pcommon);	// CSV 파일로 출력
	//<\Calculator/ShowData/FPrintf Functions>
	//<Simple calculated data : common>
void calculate_FB(void);	// 전역변수는 따로 입력값으로 쓰지 않아도 된다.
void calculate_Dvi(void);
void calculate_Dvo(void);
void calculate_Doi(void);
void calcuate_angle(void);
void calculate_e(void);
void calculate_DTh(void);
void initialize_copper(void);
void calculate_Coldtemps(ProData* prodata, Tube* ptube, Annular* pannular, int index);
	//<\Simple calculated data : common>
	//<Simple calculated data : tube>	반환값 전부다 double로
double cal_tube_Pr(double Cp, double viscosity, double conductivity);	//  모두 변수 필요.
double cal_tube_e_ratio(void); // only 전역변수
double cal_tube_p_ratio(void);// only 전역변수
double cal_tube_angle_ratio(void);// only 전역변수
double cal_tube_velocity(double density);
double cal_tube_Re(double velocity, double density, double viscosity);
double cal_tube_Nu(double Re, double e_ratio, double p_ratio, double angle_ratio, double Pr);
double cal_tube_fi(double Re, double e_ratio, double p_ratio, double angle_ratio);
double cal_tube_hi(double Nu, double conductivity);	//★가장 마지막에 구할 것.
double cal_tube_DP(double fi, double length, double density, double velocity);	//★가장 마지막에 구할 것.
	//<\Simple calculated data : tube>
	//<Simple calculated data : annular>
double cal_annular_Pr(double Cp, double viscosity, double conductivity);
double cal_annular_e_ratio(void); // only 전역변수
double cal_annular_p_ratio(void); // only 전역변수
double cal_annular_angle_ratio(void); // only 전역변수
double cal_annular_r_ratio(void); // only 전역변수
double cal_annular_Dhyd(void); // only 전역변수
double cal_annular_Aeff(void); // only 전역변수
double cal_annular_velocity(double density, double Aeff);
double cal_annular_Re(double velocity, double density, double Dhyd, double viscosity);
double cal_annular_ef(double Re, double e_ratio, double p_ratio, double angle_ratio, double r_ratio);
double cal_annular_fo(double Re, double r_ratio, double ef);
double cal_annular_Nu(double fo, double Re, double Pr, double e_ratio, double p_ratio, double r_ratio);
double cal_annular_ho(double Nu, double conductivity, double Dhyd);
double cal_annular_DP(double fo, double length, double Dhyd, double density, double velocity);
	//<\Simple calculated data : annular>
	//<Simple calculated data : objects>
double cal_common_DTln(double hot_in, double hot_out, double cold_in, double cold_out);
double cal_common_Cucond(double hot_avg, double cold_avg); // 두 온도의 평균 온도에서 conductivity를 구한다. 큰 의미 없음.
double cal_common_UA(double DTln, double HX);
double cal_common_length(double hi, double ho, double Cucond, double UA);
double cal_common_height(double length);
double cal_common_volume(double height);
double cal_common_Ai(double length);
double cal_common_Ui(double UA, double Ai);
double cal_common_Ao(double length);
double cal_common_Uo(double UA, double Ao);
	//<\Simple calculated data : objects>
	//<Evaluate Sum of Avg Values>
void SumNAvg_Common(Common* pcommon);
void SumNAvg_Tube(Tube* ptube, Common* pcommon);
void SumNAvg_Annular(Annular* pannular, Common* pcommon);
	//<\Evaluate Sum of Avg Values>
//<\Functions>
/////////////////////////////////////////////////////////////////////////////////////
//##################################################
int main(void) {
	int loop;
	char input;
	char filename[20] = "initial.csv";
	//폴더 생성 및 작업 디렉토리 설정
	char strFolderPath[100] = {"C:\\Temp\\HX_counterflow"};
	int nResult = mkdir(strFolderPath);
	
	if( nResult == 0 ) {
		printf("폴더 생성 성공. 폴더 저장된 디렉토리: %s\n", strFolderPath);
	}
	else if(nResult == -1) {
		perror("폴더 생성 실패 - 폴더가 이미 있거나 부정확함");
		printf("Error massager : %d\n", errno);
	}
	
	nResult = chdir(strFolderPath);
	if(nResult == 0) {
		printf("이동 성공. 현재 작업 디렉토리 : %s\n", strFolderPath);
	}
	else if(nResult == -1) {
		perror("이동 실패. 현재 실행 중인 exe와 같은 디렉토리에 파일이 저장이 됩니다.");
		printf("Error massager : %d\n", errno);
	}
	
	//....폴더 생성 및 작업 디렉토리 설정
	
	printf("======================================\n");
	printf("======Counter Flow HX Calculator======\n");
	printf("Process : Geometric boundary condition...\n");
	initialize_copper();
	for(loop = 0;loop<4;loop++){ //option이 0~3
		write_given_geodata(loop);
	}
	while(1) {
		calculate_FB();	// 전역변수는 따로 입력값으로 쓰지 않아도 된다.
		calculate_Dvi();
		calculate_Dvo();
		calculate_Doi();
		calcuate_angle();
		calculate_e();
		
		Print_given_geodata();
		printf("Dbi<Deo<Doi 이고 Ns가 정수인지 등등 반드시 확인을 해주시기 바랍니다.\n");
		printf("Dbi 수정하려면 0, Deo는 1, Ns는 2, Pitch는 3을 입력하세요.\n");
		printf("그냥 넘어가려면 0~3을 제외한 아무 값을 입력하세요.\n");
		getchar();
		scanf("%c", &input);
		if('0' <= input && input <= '3'){
			loop = atoi(&input);	//엉성한 부분
			write_given_geodata(loop);
		}
		else
			break;
	}
	printf("Process Complete: Geometric boundary condition!\n\n");
	
	printf("Process : Thermal boundary condition...\n");
	for(loop=0;loop<4;loop++) { // option이 0~3
		write_given_tempdata(loop);
	} 
	
	while(1) {
		calculate_DTh();
		Print_given_tempdata();
		printf("올바르게 입력했는지 확인을 해주시기 바랍니다.\n");
		printf("H_i_T을 수정하려면 0, H_o_T는 1, C_i_T는 2, 구간갯수는 3을 입력하세요\n");
		printf("그냥 넘어가시려면 0~3을 제외한 아무 값을 입력하세요.\n");
		getchar();
		scanf("%c", &input);
		if('0' <= input && input <= '3'){
			loop = atoi(&input);
			write_given_tempdata(loop);
		}
		else
			break;
	}
	printf("Process Complete : Thermal boundary condition!\n\n");
	
	printf("Process : Creating Data Table...\n");
	
	int add = 0; // Table 작성을 위해 추가 인덱스 양 계산.
	double gap = DTh/GRID; 	//(Hot_in_temp-Hot_out_temp)/(GRID*partition);	// table 쪼개는 온도 단위. 얘는 초기에 설정 함수로 옮기자.
	double temp = Hot_out_temp;
	while(temp > Cold_in_temp) {
		temp -= gap;
		add++;
	}
	
	ProData prodata_table[GRID*partition+add+1];	// partition 나눈 만큼 hot avg 값을 정확하게 하기 위해 반으로 쪼갠다.
	Tube tube_table[partition+1];
	Annular annular_table[partition+1];
	Common common_table[partition+1];
	
	printf("Process Complete: Creating Data Table!\n\n");
	
	printf("Process : Property of Water with Temperature\n");
	property_calculator(prodata_table, add, temp, gap);
	create_property_csv(prodata_table, add);
	printf("Process Complete : Property of Water with Temperature!\n\n");
	
	printf("======================================\n");
	printf("Main Process Calculating...\n");
	Calculate_HX_partly(prodata_table, tube_table, annular_table, common_table);
	SumNAvg_Common(common_table);
	SumNAvg_Tube(tube_table, common_table);
	SumNAvg_Annular(annular_table, common_table);
	create_counterflow_csv(filename, tube_table, annular_table, common_table);
	printf("\nMain Process Finish!!\n");
	printf("======================================\n\n");
	printf("\n");
	
	
	printf("======================================\n");
	printf("Process : Parametric Analysis\n");
	for(loop = 0; loop< 4; loop++) {
		write_given_parametric_factor(loop);
	}
	while(1) {
		size= (max_factor-min_factor)*100+1;
		Print_given_parametric_factor(size);
		printf("min_factor, max_factor는 반드시 소수 둘째 자리까지만 허용됩니다. 그 아래 자릿수부터는 오류가 발생합니다.\n");
		printf("min_factor를 수정하려면 0, max_factor는 1, min_Ns는 2, max_Ns는 3을 입력하세요.\n");
		printf("그냥 넘어가시려면 0~3을 제외한 아무 값이나 입력하세요.\n");
		getchar();
		scanf("%c", &input);
		if('0' <= input &&input <= '3') {
			loop = atoi(&input);
			write_given_parametric_factor(loop);
		}
		else
			break;
	}
	// 과정을 계속 반복하기 위해서는 초기값들을 저장해놓아야 한다. 
	Dbi_ini = Dbi;
	Deo_ini = Deo;
	Doi_ini = Doi;
	Pitch_ini = Pitch;
	Ns_ini = Ns;
	//이 이후에 Dbi Deo Doi Ns Pitch의 factor 조정하는 과정을 만들어 위 과정을 다시 반복하게 만들어야 한다.
	// 저장할 데이터들 호출
	Tube Tube_Modified_Dbi[size+1];
	Tube Tube_Modified_Deo[size+1];
	Tube Tube_Modified_Doi[size+1];
	Tube Tube_Modified_Pitch[size+1];
	Tube Tube_Modified_Ns[max_Ns-min_Ns+1];
	
	Annular Annular_Modified_Dbi[size+1];
	Annular Annular_Modified_Deo[size+1];
	Annular Annular_Modified_Doi[size+1];
	Annular Annular_Modified_Pitch[size+1];
	Annular Annular_Modified_Ns[max_Ns-min_Ns+1];	
	
	Common Common_Modified_Dbi[size+1];
	Common Common_Modified_Deo[size+1];
	Common Common_Modified_Doi[size+1];
	Common Common_Modified_Pitch[size+1];
	Common Common_Modified_Ns[max_Ns-min_Ns+1];
	
	char filename1[20] = "Modify_Dbi.csv";
	for(loop = 0;loop<size+1;loop++) {	// Dbi
		// 매번 geometric condition 갱신
		Dbi = (min_factor+0.01*loop)*Dbi_ini;
		calculate_FB();
		calculate_Dvi();
		calculate_Dvo();
		calculate_Doi();
		calcuate_angle();
		calculate_e();
		
		if(Dbi +2*tw< Deo) {	//이 조건이 충족되지 않으면 해당 데이터는 모두 0으로 처리된다.
			Calculate_HX_partly(prodata_table, tube_table, annular_table, common_table);
			SumNAvg_Common(common_table);
			SumNAvg_Tube(tube_table, common_table);
			SumNAvg_Annular(annular_table, common_table);
			
			Tube_Modified_Dbi[loop] = tube_table[partition];
			Annular_Modified_Dbi[loop] = annular_table[partition];
			Common_Modified_Dbi[loop] = common_table[partition];
		}
	}
	Initialize_Geoconditions();
	//filename 때문에 문서 작성할 함수를 호출해야함.
	create_counterflow_csv(filename1, Tube_Modified_Dbi, Annular_Modified_Dbi, Common_Modified_Dbi);
	
	char filename2[20] = "Modify_Deo.csv";
	for(loop = 0;loop<size+1;loop++) { // Deo
		Deo= (min_factor+0.01*loop)*Deo_ini;
		calculate_FB();
		calculate_Dvi();
		calculate_Dvo();
		//calculate_Doi(); 일단 설정을 꺼놓자.
		calcuate_angle();
		calculate_e();
		
		if(Dbi+2*tw < Deo && Deo <= Doi) {
			Calculate_HX_partly(prodata_table, tube_table, annular_table, common_table);
			SumNAvg_Common(common_table);
			SumNAvg_Tube(tube_table, common_table);
			SumNAvg_Annular(annular_table, common_table);
			
			Tube_Modified_Deo[loop] = tube_table[partition];
			Annular_Modified_Deo[loop] = annular_table[partition];
			Common_Modified_Deo[loop] = common_table[partition];
		}
	}
	Initialize_Geoconditions();	
	//filename 때문에 문서 작성할 함수를 호출해야함.
	create_counterflow_csv(filename2, Tube_Modified_Deo, Annular_Modified_Deo, Common_Modified_Deo);
	
	char filename3[20] = "Modify_Doi.csv";
	for(loop = 0;loop<size+1;loop++) {	//Doi
		Doi = (min_factor+0.01*loop)*Doi_ini;
		calculate_FB();
		calculate_Dvi();
		calculate_Dvo();
		calcuate_angle();
		calculate_e();
		
		if(Deo <= Doi) {
			Calculate_HX_partly(prodata_table, tube_table, annular_table, common_table);
			SumNAvg_Common(common_table);
			SumNAvg_Tube(tube_table, common_table);
			SumNAvg_Annular(annular_table, common_table);
			
			Tube_Modified_Doi[loop] = tube_table[partition];
			Annular_Modified_Doi[loop] = annular_table[partition];
			Common_Modified_Doi[loop] = common_table[partition];
		}
	}
	Initialize_Geoconditions();
	//filename 때문에 문서 작성할 함수를 호출해야함.
	create_counterflow_csv(filename3, Tube_Modified_Doi, Annular_Modified_Doi, Common_Modified_Doi);
	
	char filename4[20] = "Modify_Pitch.csv";
	for(loop = 0;loop<size+1;loop++) { // Pitch
		Pitch = (min_factor+0.01*loop)*Pitch_ini;
		calculate_FB();
		calculate_Dvi();
		calculate_Dvo();
		calculate_Doi();
		calcuate_angle();
		calculate_e();
		
		Calculate_HX_partly(prodata_table, tube_table, annular_table, common_table);
		SumNAvg_Common(common_table);
		SumNAvg_Tube(tube_table, common_table);
		SumNAvg_Annular(annular_table, common_table);
		
		Tube_Modified_Pitch[loop] = tube_table[partition];
		Annular_Modified_Pitch[loop] = annular_table[partition];
		Common_Modified_Pitch[loop] = common_table[partition];
	}
	Initialize_Geoconditions();
	//filename 때문에 문서 작성할 함수를 호출해야함.
	create_counterflow_csv(filename4, Tube_Modified_Pitch, Annular_Modified_Pitch, Common_Modified_Pitch);
	
	char filename5[20] = "Modify_Ns.csv";
	for(loop = min_Ns; loop < max_Ns+1;loop++) { // Ns
		Ns = loop;
		
		calculate_FB();
		calculate_Dvi();
		calculate_Dvo();
		calculate_Doi();
		calcuate_angle();
		calculate_e();
		/*printf("\n 중간 점검\n");
		printf("Ns : %d일 때,", Ns);
		Print_given_geodata();
		system("pause");*/
		
		Calculate_HX_partly(prodata_table, tube_table, annular_table, common_table);
		SumNAvg_Common(common_table);
		SumNAvg_Tube(tube_table, common_table);
		SumNAvg_Annular(annular_table, common_table);
		
		Tube_Modified_Ns[loop-1] = tube_table[partition];
		Annular_Modified_Ns[loop-1] = annular_table[partition];
		Common_Modified_Ns[loop-1] = common_table[partition];
	}
	Initialize_Geoconditions();
	//그럼 이 파일에는 최소한 Filename을 생성하는 변수가 필요하다.
	create_counterflow_csv(filename5, Tube_Modified_Ns, Annular_Modified_Ns, Common_Modified_Ns);
	printf("\nProcess Complete: Parametric Analysis\n");
	printf("유의사항\n");
	printf("QNANE값이 뜨면 DTln을 구할 때, log(0)이 되지 않았는지 확인해주세요.\n");
	printf("Modified_...csv 파일들은 0만 입력되어있는 줄이 있습니다. 그 줄은 geometric condition과 충돌하여 계산을 하지 않은 영역입니다.\n");
	system("pause");
	return 0;
}
//########################################################
///////////////////////////////////////////////////////////////////////////////////////////////
//<Function contents> 
	//<Thermal and Geometric given data>
void write_given_geodata(int option) {
	int i = option;
	switch(i) {
		case 0:
		while(1) {
			printf("Dbi 값을 입력해주십시오[m]. (ex. 15 mm => 1.5E-02) : ");
			scanf("%lf", &Dbi);
			if(Dbi < 0) {
				printf("Error : 적절하지 못한 Dbi값입니다!!\n");
			}
			else
				break;
		}	
			return;
		case 1:
			while(1) {
				printf("Deo 값을 입력해주십시오[m]. (ex. 25 mm => 2.5E-02) : ");
				scanf("%lf", &Deo);
				if(Deo < Dbi) {
					printf("Error : Deo가 Dbi보다 작으면 안됩니다!!\n");
					printf("Dbi : %lf[m], Deo : %lf[m]\n", Dbi, Deo);
				}
				else
					break;
			}			
			return;
		case 2:
			while(1) {
				printf("Ns 값을 입력해주십시오[int]. (ex. 4) : ");
				scanf("%d", &Ns);			
				if(Ns < 1) {
					printf("적절치 못한 Ns 값입니다!!\n");
				}
				else 
					break;
			}
			return;
		case 3:
			while(1) {
				printf("Pitch 값을 입력해주십시오[m]. (ex. 50/m =>2.0E-02) : ");
				scanf("%lf", &Pitch);			
				if(Pitch < 0) {
					printf("Error : 적절하지 못한 피치거리입니다!!\n");
				}
				else
					break;
			}
			return;
	}
}

void Print_given_geodata(void) {
	printf("----------------------------------------------\n");
	printf("Geometric given data 최종 확인 작업입니다.\n");
	printf("Dbi:%lf m\n", Dbi);
	printf("Deo:%lf m\n", Deo);
	printf("Ns:%d \n", Ns);
	printf("Pitch:%lf m\n", Pitch);
	printf("FB:%lf\n", FB);
	printf("Dvi:%lf\n", Dvi);
	printf("Dvo:%lf\n", Dvo);
	printf("Doi:%lf\n", Doi);
	printf("angle:%lf\n", angle);
	printf("e:%lf\n", e);
	printf("----------------------------------------------\n");
}

void write_given_tempdata(int option) {
	int i = option;
	switch(i) {
		case 0:	// hot in
			printf("Hot_in_temp을 입력해주세요[C] (ex. 70) : ");
			scanf("%lf", &Hot_in_temp);	// 수가 아닌 문자를 받았을 때 예외처리하는 방법을 찾아보자.
			return;
		case 1:		// hot out
			printf("Hot_out_temp을 입력해주세요[C] (ex. 45) : ");
			scanf("%lf", &Hot_out_temp);
			return;
		case 2:	// cold in
			while(1) {
				printf("Cold_in_temp을 입력해주세요[C] (ex.35) : ");
				scanf("%lf", &Cold_in_temp);			
				if(Cold_in_temp >= Hot_out_temp) {
					printf("적절하지 못한 Cold in temperature 값입니다.\n");
					printf("Cold_in_temp : %lf[C], Hot_out_temp : %lf[C]\n", Cold_in_temp, Hot_out_temp);
				}
				else
					break;
			}
			return;
		case 3:	// partition
			while(1) {
				printf("LMTD Method를 쓸 때 몇 개의 구간으로 나누겠습니까[int]? (최소 20이상) : ");
				scanf("%d", &partition);
				if(partition < 20 ) {
					printf("숫자가 너무 작습니다.\n");
				}
				else if(partition*GRID > 100000) {
					printf("숫자가 너무 큽니다. 구간 갯수가 %.0lf을 넘지 않도록 주의해주세요.\n", (double)100000/GRID);
				}
				else
					break;
			}
			return;
	}
}

void Print_given_tempdata(void) {
	printf("Thermal given data 최종 확인 과정입니다.\n");
	printf("----------------------------------------------\n");
	printf("Hot_in_temp : %lf\n", Hot_in_temp);
	printf("Hot_out_temp:%lf\n", Hot_out_temp);
	printf("Cold_in_temp:%lf\n", Cold_in_temp);
	printf("n개의 구간:%d\n",partition);
	printf("DTh :%lf\n", DTh);
	printf("----------------------------------------------\n");
	
}
void write_given_parametric_factor(int option) {
	int i = option;
	switch(i) {
		case 0 :
			while(1) {
				printf("Dbi, Deo, Doi, Pitch에 대한 min factor를 기입하세요 (ex. 0.8) : ");
				scanf("%lf", &min_factor);
				if(min_factor < 0 || min_factor > 1) {
					printf("Error : min_factor 범위는 0에서 1까지 입니다.\n");
				}
				else
					break;
			}
				return;
		case 1 :
			while(1) {
				printf("Dbi, Deo, Doi, Pitch에 대한 max factor를 기입하세요 (ex. 1.2) : ");
				scanf("%lf", &max_factor);
				if(max_factor-min_factor<0.01 || max_factor > 3.3) {
					printf("Error : max_factor는 %1.2lf에서 3.3까지입니다.\n", min_factor+0.01);
				}
				else 
					break;
			}
				return;
		case 2 :
			while(1) {
				printf("Ns의 최솟값을 입력해주세요 (ex. 1) : ");
				scanf("%d", &min_Ns);
				if(min_Ns < 1) {
					printf("Error : min_Ns는 적어도 1이상의 자연수여야만 합니다.\n");
				}
				else 
					break;
			}
				return;
		case 3 :
			while(1) {
				printf("Ns의 최댓값을 입력해주세요(ex. 8) : ");
				scanf("%d", &max_Ns);
				if(max_Ns <= min_Ns) {
					printf("Error : max_Ns는 min_Ns(: %d)보다 커야합니다.\n", min_Ns);
				}
				else
					break;
			}
				return;
	}
}

void Print_given_parametric_factor(int size) {
	printf("Parametric Analysis 들어가기 전 factor 최종확인 과정입니다.\n");
	printf("----------------------------------------------\n");
	printf("Dbi, Deo, Doi, Pitch의 min_factor 값 : %1.2lf\n", min_factor);
	printf("Dbi, Deo, Doi, Pitch의 max_factor 값 : %1.2lf\n", max_factor);
	printf("Ns 최솟값 : %d\n", min_Ns);
	printf("Ns 최댓값 : %d\n", max_Ns);
	printf("한 변수당 출력될 데이터 갯수 : %d\n", size);
	printf("----------------------------------------------\n");
}

void Initialize_Geoconditions(void) {
	Dbi = Dbi_ini;
	Deo = Deo_ini;
	Doi = Doi_ini;
	Ns = Ns_ini;
	Pitch = Pitch_ini;
}

	//<\Thermal and Geometric given data>
void create_property_csv(ProData* prodata, int add) {
	printf("\n	>Creating [Water_Property.csv] File\n");	//강제 형변환.
	int index;
	
	FILE *fp;
	fp = fopen("Water_Property.csv","w+");
	
	fprintf(fp, "Temp[C],Thermal Conductivity[W/m*K],Viscosity[kg/m*s],Density[kg/m^3],Cp[J/kg*K],Enthalpy[J/kg]\n");
	for(index =0;index<GRID*partition+add+1;index++) {
		fprintf(fp, "%lf,%E,%E,%E,%E,%E\n", prodata[index].temperature, prodata[index].conductivity, prodata[index].viscosity, prodata[index].density, prodata[index].Cp, prodata[index].enthalpy);		
		//printf("		>>processing : [%3.1f %c]\n", (float)(index*100)/(GRID*partition+add+1), '%');
	}
	printf("	>Complete the making Water_Property.csv file!!\n");
	fclose(fp);
}

void property_calculator(ProData* prodata, int add, double temp, double gap) {	//  
	int index, i;
	
	printf("	>계산시작 info_range : %lf ~%lf,	gap : %lf\n", Cold_in_temp, Hot_in_temp, gap);
	for(index=0; index<add+GRID*partition+1;index++) {
		prodata[index].temperature = temp;
		prodata[index].conductivity = 0;
		prodata[index].viscosity = 0;
		prodata[index].density=0;
		prodata[index].Cp=0;
		prodata[index].enthalpy=0;
		temp += gap;
		
		for(i=0;i<10;i++) {
			prodata[index].conductivity += T_cond[i]*pow(temp, i);
			prodata[index].viscosity += Viscosity[i]*pow(temp, i);
			prodata[index].density += Density[i]*pow(temp, i);
			prodata[index].Cp += 1000*Cp[i]*pow(temp, i);
			prodata[index].enthalpy += 1000*Enthalpy[i]*pow(temp, i); //둘다 KJ 단위가 기본인데, Pr 값 등 1000을 따로 곱하기는 번거로우니 여기서 그냥 1000을 곱한다.
		}
		//printf("		>>property 계산 temp : %2.3lf C[%3.1f%c] 진행됨\n", temp, (float)index*100/(add+GRID*partition+1), '%');
	}
	printf("	>property 계산 완료.\n");
}

void Calculate_HX_partly(ProData* prodata, Tube* ptube, Annular* pannular, Common* pcommon) {	//아직 점검이 안 된 함수 중 하나.
	int i=0, hot_ini_idx, cold_ini_idx;
	//초기화 과정.
	while(prodata[i].temperature < Cold_in_temp) {
		i++;
	}
	cold_ini_idx = i;
	
	while(prodata[i].temperature < Hot_out_temp) {
		i++;
	}
	hot_ini_idx = i;
	//printf("cold_ini_idx : %d,	hot_ini_idx : %d \n", cold_ini_idx, hot_ini_idx); // idx 확인. csv 파일에서 값을 찾을 때는 이 두 값에서 2씩 더하면 초기 table index 값이 나온다.
	
	ptube[0].index = hot_ini_idx;
	pannular[0].index = cold_ini_idx;
	ptube[0].in_temp = prodata[ptube[0].index].temperature;	// Hot 초기값 기입
	pannular[0].in_temp = prodata[pannular[0].index].temperature;	// Cold 초기값 기입
	ptube[0].in_enthalpy = prodata[ptube[0].index].enthalpy;	// Hot 사이드 초기 엔탈피
	pannular[0].in_enthalpy = prodata[pannular[0].index].enthalpy;// Cold 초기 엔탈피

	for(i=0;i<partition;i++) {	// 이 아래부터는 계속 반복과정. 좀 더 메모리 적게 쓸 수 있으나 너무 복잡해서 구별하기 쉽게 구조체에 데이터 조금 더 쓴다.
		// STEP 1: hot side Tin Tout 구하기		
		ptube[i+1].index = ptube[i].index + GRID; // Hot (i+1) 테이블 인덱스 값 입력. 얘는 간격이 일정하다. 얘들은 나중에 그 따로 마지막에서 메모리 할당 안 되서 오류 일으킬 건데 그 부분은 다시 따로 손 봐야 한다.
		ptube[i].out_temp = ptube[i].in_temp+DTh;// Hot(i) out 온도 기입.
		ptube[i+1].in_temp = ptube[i].out_temp; // i번째 out이 i+1번째 in과 동일함.
		
		//STEP2 : 열교환량 계산하기
		ptube[i].out_enthalpy = prodata[ptube[i+1].index].enthalpy; // hot(i) out 엔탈피값
		ptube[i+1].in_enthalpy = ptube[i].out_enthalpy;	// hot(i+1) in 엔탈피값 
		pcommon[i].HX = mh*(ptube[i].out_enthalpy-ptube[i].in_enthalpy);	// part i에서의 열교환량 구함.
		
		//STEP3 : Cold_out 탐색후 i+1까지 기입.
		pannular[i+1].index = pannular[i].index+1; // pannular i+1 인덱스 초기화
		while(pcommon[i].HX > mc*((prodata[pannular[i+1].index].enthalpy)-(prodata[pannular[i].index].enthalpy))) {
			pannular[i+1].index++;
		}// 자동으로 cold(i+1) 테이블 인덱스 기입.
		pannular[i].out_enthalpy = prodata[pannular[i+1].index].enthalpy;
		pannular[i+1].in_enthalpy = pannular[i].out_enthalpy;
		pannular[i].out_temp = prodata[pannular[i+1].index].temperature;
		pannular[i+1].in_temp = pannular[i].out_temp;
		
		//STEP4: 나머지 온도 관련 프로퍼티 전부 구하기
		ptube[i].avg_temp=(prodata[ptube[i].index].temperature + prodata[ptube[i+1].index].temperature)/2;
		ptube[i].conductivity=(prodata[ptube[i].index].conductivity + prodata[ptube[i+1].index].conductivity)/2;
		ptube[i].viscosity=(prodata[ptube[i].index].viscosity + prodata[ptube[i+1].index].viscosity)/2;
		ptube[i].density=(prodata[ptube[i].index].density + prodata[ptube[i+1].index].density)/2;
		ptube[i].Cp=(prodata[ptube[i].index].Cp + prodata[ptube[i+1].index].Cp)/2;
		
		pannular[i].avg_temp=(prodata[pannular[i].index].temperature + prodata[pannular[i+1].index].temperature)/2;
		pannular[i].conductivity=(prodata[pannular[i].index].conductivity + prodata[pannular[i+1].index].conductivity)/2;
		pannular[i].viscosity=(prodata[pannular[i].index].viscosity + prodata[pannular[i+1].index].viscosity)/2;
		pannular[i].density=(prodata[pannular[i].index].density + prodata[pannular[i+1].index].density)/2;
		pannular[i].Cp=(prodata[pannular[i].index].Cp + prodata[pannular[i+1].index].Cp)/2;		
		//위 코드들 점검용 출력 코드.
		/*printf("hot[table idx : %d] : %lf C, %lf W/mK, %lf Pa*s, %lf kg/m^3, %lf kJ/kg*K\n", ptube[i].index, ptube[i].avg_temp, ptube[i].conductivity, ptube[i].viscosity, ptube[i].density, ptube[i].Cp);
		printf("cold[table idx : %d] : %lf C, %lf W/mK, %lf Pa*s, %lf kg/m^3, %lf kJ/kg*K\n", pannular[i].index, pannular[i].avg_temp, pannular[i].conductivity, pannular[i].viscosity, pannular[i].density, pannular[i].Cp);
		printf("HX[%d 번째] : %lf \n", i, pcommon[i].HX);*/
		
		//위 까지는 반드시! 구조체의 온도 관련 프로퍼티들을 전부 채워야 한다. in_temp, index, in_enthalpy는 i+1관련 기입을 했는지 확인할 것.
		// step4 까지는 이상 무. csv 출력 파일이나 만들자.
		//STEP 5 : tube side 주요 값들 구하기 hi 값과 DP은 가장 마지막에 구한다.
		ptube[i].Pr = cal_tube_Pr(ptube[i].Cp, ptube[i].viscosity, ptube[i].conductivity);
		ptube[i].e_ratio = cal_tube_e_ratio();
		ptube[i].p_ratio = cal_tube_p_ratio();
		ptube[i].angle_ratio = cal_tube_angle_ratio();
		ptube[i].velocity = cal_tube_velocity(ptube[i].density);
		ptube[i].Re = cal_tube_Re(ptube[i].velocity, ptube[i].density, ptube[i].viscosity);
		ptube[i].Nu = cal_tube_Nu(ptube[i].Re, ptube[i].e_ratio, ptube[i].p_ratio, ptube[i].angle_ratio, ptube[i].Pr);
		ptube[i].hi = cal_tube_hi(ptube[i].Nu, ptube[i].conductivity);
		ptube[i].fi = cal_tube_fi(ptube[i].Re, ptube[i].e_ratio, ptube[i].p_ratio, ptube[i].angle_ratio);
		
		//STEP 6: annular side의 주요 값들 구하기 ho 값과 DP는 나중에 구한다.
		pannular[i].Pr = cal_annular_Pr(pannular[i].Cp, pannular[i].viscosity, pannular[i].conductivity);
		pannular[i].e_ratio = cal_annular_e_ratio();
		pannular[i].p_ratio = cal_annular_p_ratio();
		pannular[i].angle_ratio = cal_annular_angle_ratio();
		pannular[i].r_ratio = cal_annular_r_ratio();
		pannular[i].Dhyd = cal_annular_Dhyd();
		pannular[i].Aeff = cal_annular_Aeff();
		pannular[i].velocity = cal_annular_velocity(pannular[i].density, pannular[i].Aeff);
		pannular[i].Re = cal_annular_Re(pannular[i].velocity, pannular[i].density, pannular[i].Dhyd, pannular[i].viscosity);
		pannular[i].ef = cal_annular_ef(pannular[i].Re, pannular[i].e_ratio, pannular[i].p_ratio, pannular[i].angle_ratio, pannular[i].r_ratio);
		pannular[i].fo = cal_annular_fo(pannular[i].Re, pannular[i].r_ratio, pannular[i].ef);
		pannular[i].Nu = cal_annular_Nu(pannular[i].fo, pannular[i].Re, pannular[i].Pr, pannular[i].e_ratio, pannular[i].p_ratio, pannular[i].r_ratio);
		pannular[i].ho = cal_annular_ho(pannular[i].Nu, pannular[i].conductivity, pannular[i].Dhyd);
		
		//STEP 7: UA, HX, Length 등 common side의 값들을 전부 구한다. 
		//pcommon[i].DTln = ((ptube[i].in_temp-pannular[i].in_temp)-(ptube[i].out_temp-ptube[i].out_temp))/ln((ptube[i].in_temp-pannular[i].in_temp)/(ptube[i].out_temp-ptube[i].out_temp));  의미 없는 개노가다..
		pcommon[i].DTln = cal_common_DTln(ptube[i].in_temp, ptube[i].out_temp, pannular[i].in_temp, pannular[i].out_temp);
		pcommon[i].Cucond = cal_common_Cucond(ptube[i].avg_temp, pannular[i].avg_temp);
		pcommon[i].UA = cal_common_UA(pcommon[i].DTln, pcommon[i].HX);
		pcommon[i].length = cal_common_length(ptube[i].hi, pannular[i].ho, pcommon[i].Cucond, pcommon[i].UA);
		pcommon[i].height = cal_common_height(pcommon[i].length);
		pcommon[i].volume = cal_common_volume(pcommon[i].height);
		pcommon[i].Ai = cal_common_Ai(pcommon[i].length);
		pcommon[i].Ui = cal_common_Ui(pcommon[i].UA, pcommon[i].Ai);
		pcommon[i].Ao = cal_common_Ao(pcommon[i].length);
		pcommon[i].Uo = cal_common_Uo(pcommon[i].UA, pcommon[i].Ao);
		
		//STEP 8 : DPi DPo 값 구하기
		ptube[i].DP = cal_tube_DP(ptube[i].fi, pcommon[i].length, ptube[i].density, ptube[i].velocity);
		pannular[i].DP = cal_annular_DP(pannular[i].fo, pcommon[i].length, pannular[i].Dhyd, pannular[i].density, pannular[i].velocity);
		}
}

void create_counterflow_csv(char* filename, Tube* ptube, Annular* pannular, Common* pcommon) {	// Filename도 받아서 적절히 이름도 바꿀 수 있게 만들자.
	int index;
	char name[80] = "Conditions_";
	char name1[80] = "CounterFlow_LMTD_Tubeside_";
	char name2[80] = "CounterFlow_LMTD_Annularside_";
	char name3[80] = "CounterFlow_LMTD_Results_";
	FILE* fp;
	
	if(!strcmp("initial.csv", filename)) {
		strcat(name, filename);
		printf("\n	>Creating[%s] File\n", name);
		fp = fopen(name, "w+");
		
		fprintf(fp, "Geometric Conditions\n");
		fprintf(fp, "Dbi[m],Deo[m],tw[m],Ns,Pitch[m],Dcan[m],R,FB,Dvi[m],Dvo[m],Doi[m],angle,e\n");
		fprintf(fp, "%E,%E,%E,%d,%E,%E,%E,%E,%E,%E,%E,%E,%E\n", Dbi, Deo, tw, Ns, Pitch, Dcan, R, FB, Dvi, Dvo, Doi, angle, e);
		
		fprintf(fp, "Thermal and Fluid Conditions\n");
		fprintf(fp, "Tube Side,,,Annular Side,,,\n");
		fprintf(fp, "Hot_in_temp[C],Hot_out_temp[C],mh[kg/s],Cold_in_temp[C],Cold_out_temp[C], mc[kg/s]\n");
		fprintf(fp, "%lf, %lf, %lf, %lf, %lf, %lf\n",Hot_in_temp,Hot_out_temp,mh,Cold_in_temp,Cold_out_temp, mc);
		
		fclose(fp);
		
		//Step2 : Tubeside 데이터 csv 만들기
		double x_position = pcommon[partition].length;
		printf("\n	>Creating [CounterFlow_LMTD_Tubeside_initial] File\n");
		strcat(name1, filename);
		fp = fopen(name1, "w+");
		fprintf(fp, "index[#],in_temp[C],out_temp[C],avg_temp[C],in_enthalpy[J],out_enthalpy[J],conductivity[W/mK],visocosity[],density[kg/m^3],Cp[J/kg*K],Pr,e_ratio,p_ratio,angle_ratio,velocity,Re,Nu,hi,fi,DP,,Length[m], x_position[m]\n");
		fprintf(fp, "average & Sum values\n");
		fprintf(fp, "%d,%lf,%lf,%lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", partition,ptube[partition].in_temp,ptube[partition].out_temp,ptube[partition].avg_temp,ptube[partition].in_enthalpy,ptube[partition].out_enthalpy,ptube[partition].conductivity,ptube[partition].viscosity,ptube[partition].density,ptube[partition].Cp,ptube[partition].Pr,ptube[partition].e_ratio,ptube[partition].p_ratio,ptube[partition].angle_ratio,ptube[partition].velocity,ptube[partition].Re,ptube[partition].Nu,ptube[partition].hi,ptube[partition].fi,ptube[partition].DP,pcommon[partition].length);
		fprintf(fp, "Details!\n");
		for(index=0;index<partition;index++) {
			fprintf(fp, "%d,%lf,%lf,%lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E,%1.5E\n", index,ptube[index].in_temp,ptube[index].out_temp,ptube[index].avg_temp,ptube[index].in_enthalpy,ptube[index].out_enthalpy,ptube[index].conductivity,ptube[index].viscosity,ptube[index].density,ptube[index].Cp,ptube[index].Pr,ptube[index].e_ratio,ptube[index].p_ratio,ptube[index].angle_ratio,ptube[index].velocity,ptube[index].Re,ptube[index].Nu,ptube[index].hi,ptube[index].fi,ptube[index].DP,pcommon[index].length,x_position);
			x_position -= pcommon[index].length;
		}

		printf("\n	Complete : Creating[CounterFlow_LMTD_Tubeside_initial]File\n");
		fclose(fp);
		
		//Step3 : Annular side 데이터 csv 만들기
		x_position = pcommon[partition].length;
		printf("\n	>Creating[CounterFlow_LMTD_Annularside_initial] File\n");
		strcat(name2, filename);
		fp = fopen(name2, "w+");
		
		fprintf(fp, "index[#],in_temp[C],out_temp[C],avg_temp[C],in_enthalpy[J],out_enthalpy[J],conductivity[W/mK],visocosity[],density[kg/m^3],Cp[J/kg*K],Pr,e_ratio,p_ratio,angle_ratio,r_ratio,Dhyd,Aeff,velocity,Re,Nu,ef,ho,fo,DP,,Length[m],x_position[m]\n"); //15,17,18,21
		fprintf(fp, "average & Sum Values\n");
		fprintf(fp, "%d,%lf,%lf,%lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", partition,pannular[partition].in_temp,pannular[partition].out_temp,pannular[partition].avg_temp,pannular[partition].in_enthalpy,pannular[partition].out_enthalpy,pannular[partition].conductivity,pannular[partition].viscosity,pannular[partition].density,pannular[partition].Cp,pannular[partition].Pr,pannular[partition].e_ratio,pannular[partition].p_ratio,pannular[partition].angle_ratio,pannular[partition].r_ratio,pannular[partition].Dhyd,pannular[partition].Aeff,pannular[partition].velocity,pannular[partition].Re,pannular[partition].Nu,pannular[partition].ef,pannular[partition].ho,pannular[partition].fo,pannular[partition].DP,pcommon[partition].length);
		fprintf(fp, "Details!\n");
		for(index=0;index<partition;index++) {
			fprintf(fp, "%d,%lf,%lf,%lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E,%1.5E\n", index,pannular[index].in_temp,pannular[index].out_temp,pannular[index].avg_temp,pannular[index].in_enthalpy,pannular[index].out_enthalpy,pannular[index].conductivity,pannular[index].viscosity,pannular[index].density,pannular[index].Cp,pannular[index].Pr,pannular[index].e_ratio,pannular[index].p_ratio,pannular[index].angle_ratio,pannular[index].r_ratio,pannular[index].Dhyd,pannular[index].Aeff,pannular[index].velocity,pannular[index].Re,pannular[index].Nu,pannular[index].ef,pannular[index].ho,pannular[index].fo,pannular[index].DP,pcommon[index].length,x_position);
			x_position -= pcommon[index].length;
		}
		printf("\n	Complete : Creating[CounterFlow_LMTD_Annularside_initial]File\n");
		fclose(fp);
		
		//Step4 : Common side 데이터 csv 만들기
		printf("\n	Creating[CounterFlow_LMTD_Results_initial] File\n");
		strcat(name3, filename);
		fp = fopen(name3, "w+");
		fprintf(fp, "index[#],HX[W],length[m],height[m],volume[m^3],DTln[K],UA[W/K],Ui[W/K*m^2],Ai[m^2],Uo[W/K*m^2],Ao[m^2],Cucond[W/mK],,length[m]\n");
		fprintf(fp, "average & Sum Values\n");
		fprintf(fp, "%d,%1.5E,%1.5E,%1.5E,%1.5E,%2.4lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E\n", partition,pcommon[partition].HX,pcommon[partition].length,pcommon[partition].height,pcommon[partition].volume,pcommon[partition].DTln,pcommon[partition].UA,pcommon[partition].Ui,pcommon[partition].Ai,pcommon[partition].Uo,pcommon[partition].Ao,pcommon[partition].Cucond);
		fprintf(fp, "Details!!\n");
		for(index = 0;index<partition;index++){
			fprintf(fp, "%d,%1.5E,%1.5E,%1.5E,%1.5E,%2.4lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E\n", index,pcommon[index].HX,pcommon[index].length,pcommon[index].height,pcommon[index].volume,pcommon[index].DTln,pcommon[index].UA,pcommon[index].Ui,pcommon[index].Ai,pcommon[index].Uo,pcommon[index].Ao,pcommon[index].Cucond);
		}
		printf("\n	Complete : Creating[CounterFlow_LMTD_Results_initial]File\n");
		fclose(fp);
		// 추가로 최종 평균 데이터를 담고 있는 csv파일을 저장해야 한다.
	}
	else if (!strcmp("Modify_Ns.csv", filename)) {	//Ns는 따로 해줘야 한다. 범위가 달라서.
		strcat(name, filename);
		printf("\n	>Creating[%s] File\n", name);
		fp = fopen(name, "w+");
		
		/*fprintf(fp, "Geometric Conditions\n");
		fprintf(fp, "Dbi[m],Deo[m],tw[m],Ns,Pitch[m],Dcan[m],R,FB,Dvi[m],Dvo[m],Doi[m],angle,e\n");
		fprintf(fp, "%E,%E,%E,%d,%E,%E,%E,%E,%E,%E,%E,%E,%E\n", Dbi, Deo, tw, Ns, Pitch, Dcan, R, FB, Dvi, Dvo, Doi, angle, e);
		*/	 // 시간 부족으로 geometric condition이 팩터별로 변하는 걸 구현 못함.
		fprintf(fp, "Thermal and Fluid Conditions\n");
		fprintf(fp, "Tube Side,,,Annular Side,,,\n");
		fprintf(fp, "Hot_in_temp[C],Hot_out_temp[C],mh[kg/s],Cold_in_temp[C],Cold_out_temp[C], mc[kg/s]\n");
		fprintf(fp, "%lf, %lf, %lf, %lf, %lf, %lf\n",Hot_in_temp,Hot_out_temp,mh,Cold_in_temp,Cold_out_temp, mc);
		
		fclose(fp);
		/////////
		strcat(name1, filename);
		printf("\n	>Creating [%s] File\n", name1);
		fp= fopen(name1, "w+");
		
		fprintf(fp, "index[#],conductivity[W/mK],viscosity[kg/m*s],density[kg/m^3],Cp[J/kg*K],Pr,e_ratio,p_ratio,angle_ratio,velocity,Re,Nu,hi,fi,DP,,length[m]\n");
		for(index = 0;index<max_Ns-min_Ns+1;index++) {
			fprintf(fp, "%d,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", min_Ns+index,ptube[index].conductivity,ptube[index].viscosity,ptube[index].density,ptube[index].Cp,ptube[index].Pr,ptube[index].e_ratio,ptube[index].p_ratio,ptube[index].angle_ratio,ptube[index].velocity,ptube[index].Re,ptube[index].Nu,ptube[index].hi,ptube[index].fi,ptube[index].DP,pcommon[index].length);
		}
		printf("\n	>Complete : Creating[%s] File\n", name1);
		fclose(fp);
		//////////////////
		strcat(name2,filename);
		printf("\n	>Creating [%s] File\n", name2);
		fp = fopen(name2, "w+");
	
		fprintf(fp, "index[#],conductivity[W/mK],viscosity[kg/m*s],density[kg/m^3],Cp[J/kg*K],Pr,e_ratio,p_ratio,angle_ratio,r_ratio,Dhyd,Aeff,velocity,Re,Nu,ef,ho,fo,DP,,length[m]\n");
		
		for(index = 0;index<max_Ns-min_Ns+1;index++) {
			fprintf(fp, "%d,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", min_Ns+index,pannular[index].conductivity,pannular[index].viscosity,pannular[index].density,pannular[index].Cp,pannular[index].Pr,pannular[index].e_ratio,pannular[index].p_ratio,pannular[index].angle_ratio,pannular[index].r_ratio,pannular[index].Dhyd,pannular[index].Aeff,pannular[index].velocity,pannular[index].Re,pannular[index].Nu,pannular[index].ef,pannular[index].ho,pannular[index].fo,pannular[index].DP,pcommon[index].length);
		}
		printf("\n	Complete : Creating[%s] File\n", name2);
		fclose(fp);
		/////////////////////
		strcat(name3, filename);
		printf("\n	>Creating [%s] File\n", name3);
		fp=fopen(name3, "w+");

		fprintf(fp, "index[#],HX[W],length[m],height[m],volume[m^3],DTln[K],UA[W/K],Ui[W/K*m^2],Ai[m^2],Uo[W/K*m^2],Ao[m^2],Cucond[W/mK],,length[m]\n");
		
		for(index = 0; index<max_Ns-min_Ns+1;index++) {
			fprintf(fp, "%d,%1.5E,%1.5E,%1.5E,%1.5E,%2.4lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", min_Ns+index,pcommon[index].HX,pcommon[index].length,pcommon[index].height,pcommon[index].volume,pcommon[index].DTln,pcommon[index].UA,pcommon[index].Ui,pcommon[index].Ai,pcommon[index].Uo,pcommon[index].Ao,pcommon[index].Cucond,pcommon[index].length);
		}
		printf("\n	Complete : Creating [%s] File\n", name3);
		fclose(fp);
	}
	
	else {	// parametric analysis 일 경우 그냥 한 문서에 factor에 관한 최종값들만 모아놓는다.
		strcat(name, filename);
		printf("\n	>Creating[%s] File\n", name);
		fp = fopen(name, "w+");
		
		/*fprintf(fp, "Geometric Conditions\n");
		fprintf(fp, "Dbi[m],Deo[m],tw[m],Ns,Pitch[m],Dcan[m],R,FB,Dvi[m],Dvo[m],Doi[m],angle,e\n");
		fprintf(fp, "%E,%E,%E,%d,%E,%E,%E,%E,%E,%E,%E,%E,%E\n", Dbi, Deo, tw, Ns, Pitch, Dcan, R, FB, Dvi, Dvo, Doi, angle, e);*/
		
		fprintf(fp, "Thermal and Fluid Conditions\n");
		fprintf(fp, "Tube Side,,,Annular Side,,,\n");
		fprintf(fp, "Hot_in_temp[C],Hot_out_temp[C],mh[kg/s],Cold_in_temp[C],Cold_out_temp[C], mc[kg/s]\n");
		fprintf(fp, "%lf, %lf, %lf, %lf, %lf, %lf\n",Hot_in_temp,Hot_out_temp,mh,Cold_in_temp,Cold_out_temp, mc);
		
		fclose(fp);
		////////
		strcat(name1, filename);
		printf("\n	>Creating [%s] File\n", name1);
		fp= fopen(name1, "w+");
		fprintf(fp, "index[#],conductivity[W/mK],viscosity[kg/m*s],density[kg/m^3],Cp[J/kg*K],Pr,e_ratio,p_ratio,angle_ratio,velocity,Re,Nu,hi,fi,DP,, length[m]\n");
		for(index = 0;index<size+1;index++) {
			fprintf(fp, "%lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", min_factor+0.01*index,ptube[index].conductivity,ptube[index].viscosity,ptube[index].density,ptube[index].Cp,ptube[index].Pr,ptube[index].e_ratio,ptube[index].p_ratio,ptube[index].angle_ratio,ptube[index].velocity,ptube[index].Re,ptube[index].Nu,ptube[index].hi,ptube[index].fi,ptube[index].DP,pcommon[index].length);
		}
		printf("\n	>Complete : Creating[%s] File\n", name1);
		fclose(fp);
		//////////////////
		strcat(name2,filename);
		printf("\n	>Creating [%s] File\n", name2);
		fp = fopen(name2, "w+");
		
		fprintf(fp, "index[#],conductivity[W/mK],viscosity[kg/m*s],density[kg/m^3],Cp[J/kg*K],Pr,e_ratio,p_ratio,angle_ratio,r_ratio,Dhyd,Aeff,velocity,Re,Nu,ef,ho,fo,DP,,length[m]\n");
		
		for(index = 0;index<size+1;index++) {
			fprintf(fp, "%lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", min_factor+0.01*index,pannular[index].conductivity,pannular[index].viscosity,pannular[index].density,pannular[index].Cp,pannular[index].Pr,pannular[index].e_ratio,pannular[index].p_ratio,pannular[index].angle_ratio,pannular[index].r_ratio,pannular[index].Dhyd,pannular[index].Aeff,pannular[index].velocity,pannular[index].Re,pannular[index].Nu,pannular[index].ef,pannular[index].ho,pannular[index].fo,pannular[index].DP,pcommon[index].length);
		}
		printf("\n	Complete : Creating[%s] File\n", name2);
		fclose(fp);
		/////////////////////
		strcat(name3, filename);
		printf("\n	>Creating [%s] File\n", name3);
		fp=fopen(name3, "w+");

		fprintf(fp, "index[#],HX[W],length[m],height[m],volume[m^3],DTln[K],UA[W/K],Ui[W/K*m^2],Ai[m^2],Uo[W/K*m^2],Ao[m^2],Cucond[W/mK],,length[m]\n");
		
		for(index = 0; index<size+1;index++) {
			fprintf(fp, "%lf,%1.5E,%1.5E,%1.5E,%1.5E,%2.4lf,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,%1.5E,,%1.5E\n", min_factor+0.01*index,pcommon[index].HX,pcommon[index].length,pcommon[index].height,pcommon[index].volume,pcommon[index].DTln,pcommon[index].UA,pcommon[index].Ui,pcommon[index].Ai,pcommon[index].Uo,pcommon[index].Ao,pcommon[index].Cucond,pcommon[index].length);
		}
		printf("\n	Complete : Creating [%s] File\n", name3);
		fclose(fp);
	}
} 

	//<Simple calculated data : common>
void calculate_FB(void){
	FB = (1-R)*M_PI*Dbi/Ns;
}	// 전역변수는 따로 입력값으로 쓰지 않아도 된다.
void calculate_Dvi(void){
	Dvi = sqrt(pow(Dbi,2.0)+(Ns*(Deo-Dbi-2*tw)*FB)/M_PI);
}
void calculate_Dvo(void){
	Dvo = Dvi+2*tw;
}
void calculate_Doi(void){
	Doi = Deo+1.5E-04;
}
void calcuate_angle(void){
	angle = atan((M_PI*Dvo)/(Ns*Pitch))*180/M_PI;
}
void calculate_e(void){
	e=(Deo-(Dbi+2*tw))/2;
}
void calculate_DTh(void) {
	//double sub = Hot_in_temp-Hot_out_temp;
	DTh = (Hot_in_temp - Hot_out_temp)/partition; 	// 1.#INF00해결. partition 형 double로 해놓고 %d로 값받아서 나온 현상.
}
void initialize_copper(void) {
	copper[0].temperature = -73;
	copper[0].conductivity = 413;
	copper[1].temperature = 0;
	copper[1].conductivity = 401;
	copper[2].temperature = 127;
	copper[2].conductivity = 392;
	copper[3].temperature = 327;
	copper[3].conductivity = 383;
	copper[4].temperature = 527;
	copper[4].conductivity = 371;
	copper[5].temperature = 727;
	copper[5].conductivity = 357;
	copper[6].temperature = 927;
	copper[6].conductivity = 342;
}
	//<\Simple calculated data : common>
	//<Simple calculated data : tube>	반환값 전부다 double로
double cal_tube_Pr(double Cp, double viscosity, double conductivity){
	return Cp*viscosity/conductivity;
}
double cal_tube_e_ratio(void){
	return e/Dvi;
} 
double cal_tube_p_ratio(void){
	return Pitch/Dvi;
}
double cal_tube_angle_ratio(void){
	return angle/90;
}
double cal_tube_velocity(double density){
	return 4*mh/(M_PI*density*pow(Dvi,2.0));
}
double cal_tube_Re(double velocity, double density, double viscosity){
	return velocity*density*Dvi/viscosity;
}
double cal_tube_Nu(double Re, double e_ratio, double p_ratio, double angle_ratio, double Pr){
	return 0.064*pow(Re, 0.773)*pow(e_ratio, -0.242)*pow(p_ratio, -0.108)*pow(angle_ratio, -0.599)*pow(Pr, 0.4);
}
double cal_tube_fi(double Re, double e_ratio, double p_ratio, double angle_ratio){
	return 1.209*pow(Re, -0.261)*pow(e_ratio,1.26-0.05*p_ratio)*pow(p_ratio, -1.660+2.033*e_ratio)*pow(angle_ratio, -2.699+3.67*e_ratio);
}
double cal_tube_hi(double Nu, double conductivity){
	return Nu*conductivity/Dvi;
}
double cal_tube_DP(double fi, double length, double density, double velocity){
	return fi*length/(2*Dvi)*density*pow(velocity, 2.0);
}
	//<\Simple calculated data : tube>
	//<Simple calculated data : annular>
double cal_annular_Pr(double Cp, double viscosity, double conductivity) {
	return Cp*viscosity/conductivity;
}
double cal_annular_e_ratio(void) {
	return e/Dvo;
} // only 전역변수
double cal_annular_p_ratio(void){
	return Pitch/Dvo;
} // only 전역변수
double cal_annular_angle_ratio(void){
	return angle/90;
}
double cal_annular_r_ratio(void){
	return Dvo/Doi;
} // only 전역변수
double cal_annular_Dhyd(void){
	return Doi-Dvo;
} // only 전역변수
double cal_annular_Aeff(void){
	return M_PI/4*(pow(Doi, 2.0)-pow(Dvo, 2.0));
} // only 전역변수
double cal_annular_velocity(double density, double Aeff){
	return mc/(density*Aeff);
}
double cal_annular_Re(double velocity, double density, double Dhyd, double viscosity){
	return velocity*density*Dhyd/viscosity;
}
double cal_annular_ef(double Re, double e_ratio, double p_ratio, double angle_ratio, double r_ratio){
	return (1+222*pow(Re,0.09)*pow(e_ratio,2.4)*pow(p_ratio,-0.49)*pow(angle_ratio,-0.38)*pow(r_ratio,2.22));
}
double cal_annular_fo(double Re, double r_ratio, double ef) {
	return 4*pow(1.7372*log(Re/(1.964*log(Re)-3.8215)),-2)*(1+0.0925*r_ratio)*ef;
}
double cal_annular_Nu(double fo, double Re, double Pr, double e_ratio, double p_ratio, double r_ratio) {	// v 0.9 버전에서 발생한 오류 주 원인.
	//A/B*C구조
	double fo_modi = fo/8;
	double A_part = fo_modi*Re*Pr;
	double B_part = (1+9.77*pow(fo_modi, 0.5)*(pow(Pr, 0.66666666666666666)-1));	//  pow 안의 변수들의 계산은 자동으로 처리되지 않는다.
	double C_part = pow(Re,-0.20)*pow(e_ratio,-0.32)*pow(p_ratio,-0.28)*pow(r_ratio,-1.64);
	/*printf("현재 fo 값 : %E, Re 값 : %E, Pr 값 : %E, e* = %E, p* = %E, r* = %E 값입니다.\n", fo, Re, Pr, e_ratio, p_ratio, r_ratio);
	printf("\nA_part : %lf\nB_part : %lf\nC_part : %lf\nNu : %lf입니다.\n", A_part, B_part, C_part,A_part/B_part*C_part);
	system("pause");*/
	return A_part/B_part*C_part;
	//return (((fo/8)*Re*Pr)/(1+(9.77*sqrt(fo/8)*(pow(Pr,2/3)-1))))*pow(Re,-0.2)*pow(e_ratio,-0.32)*pow(p_ratio,-0.28)*pow(r_ratio,-1.64);
}
double cal_annular_ho(double Nu, double conductivity, double Dhyd){
	return Nu*conductivity/Dhyd;
}
double cal_annular_DP(double fo, double length, double Dhyd, double density, double velocity){
	return fo*length/(2*Dhyd)*density*pow(velocity, 2.0);
}
	//<\Simple calculated data : annular>
	//<Simple calculated data : objects>
double cal_common_DTln(double hot_in, double hot_out, double cold_in, double cold_out){
	return (((hot_in-cold_in)-(hot_out-cold_out))/log((hot_in-cold_in)/(hot_out-cold_out)));
}
double cal_common_Cucond(double hot_avg, double cold_avg){
	double avg_temp = (hot_avg+cold_avg)/2;
	double ratio;
	int i=0;
	while(avg_temp > copper[i].temperature) {	// 927도보다 높으면 오작동.
		i++;
	}	// 높은 인덱스 자동 반환.
	if(i==0) {
		return copper[0].temperature;
	}
	else {
		ratio = (avg_temp - copper[i-1].temperature)/(copper[i].temperature - copper[i-1].temperature);
		return copper[i-1].conductivity+ratio*(copper[i].conductivity - copper[i-1].conductivity);
	}
};
double cal_common_UA(double DTln, double HX){
	return HX/DTln;
}
double cal_common_length(double hi, double ho, double Cucond, double UA){
	return UA*((1/(hi*M_PI*Dvi))+((log(Dvo/Dvi))/(2*M_PI*Cucond))+(1/(ho*M_PI*Dvo)));
}
double cal_common_height(double length){
	return length/(M_PI*Dcan)*(Doi+2*tw);
}
double cal_common_volume(double height){ // 공식 미확인 함수
	return (Dcan+Doi+2*tw)*(Dcan+Doi+2*tw)*height*M_PI/4;
}
double cal_common_Ai(double length){
	return M_PI*Dvi*length;
}
double cal_common_Ui(double UA, double Ai){
	return UA/Ai;
}
double cal_common_Ao(double length){
	return M_PI*Dvo*length;
}
double cal_common_Uo(double UA, double Ao){
	return UA/Ao;
}
	//<\Simple calculated data : objects>
		//<Evaluate Sum of Avg Values>
void SumNAvg_Common(Common* pcommon)  {	// partition은 전역변수임. 
	int i;
	// 초기화 과정
	pcommon[partition].HX = 0;	// Sum
	pcommon[partition].length = 0;	// Sum //가장 중요한 변수.
	pcommon[partition].height = 0; // Sum
	pcommon[partition].volume = 0; // Sum
	pcommon[partition].DTln = 0;	// avg
	pcommon[partition].UA = 0; // Sum
	pcommon[partition].Ui = 0; // avg
	pcommon[partition].Ai = 0; // Sum
	pcommon[partition].Uo = 0;// avg
	pcommon[partition].Ao = 0; // Sum
	pcommon[partition].Cucond = 0; // avg
	
	for(i=0;i<partition;i++) {
	pcommon[partition].HX += pcommon[i].HX;	
	pcommon[partition].length += pcommon[i].length;	 
	pcommon[partition].height += pcommon[i].height; 
	pcommon[partition].volume += pcommon[i].volume; 
	pcommon[partition].DTln += pcommon[i].DTln * pcommon[i].length;	
	pcommon[partition].UA += pcommon[i].UA; 
	pcommon[partition].Ui += pcommon[i].Ui * pcommon[i].length; 
	pcommon[partition].Ai += pcommon[i].Ai; 
	pcommon[partition].Uo += pcommon[i].Uo * pcommon[i].length;
	pcommon[partition].Ao += pcommon[i].Ao; 
	pcommon[partition].Cucond += pcommon[i].Cucond * pcommon[i].length; 		
	}
	pcommon[partition].DTln /= pcommon[partition].length;
	pcommon[partition].Ui /= pcommon[partition].length;
	pcommon[partition].Uo /= pcommon[partition].length;
	pcommon[partition].Cucond /= pcommon[partition].length;
}

void SumNAvg_Tube(Tube* ptube, Common* pcommon) {
	int i;
	ptube[partition].in_temp = 0;	// 무의미한 값들.
	ptube[partition].out_temp = 0;
	ptube[partition].avg_temp = 0;
	ptube[partition].index = 0;
	ptube[partition].in_enthalpy = 0;
	ptube[partition].out_enthalpy = 0;
	
	ptube[partition].conductivity= 0;	//conductivity 부터 시작.
	ptube[partition].viscosity = 0;
	ptube[partition].density = 0;
	ptube[partition].Cp = 0;
	ptube[partition].Pr = 0;
	ptube[partition].e_ratio = 0;
	ptube[partition].p_ratio = 0;
	ptube[partition].angle_ratio = 0;
	ptube[partition].velocity = 0;
	ptube[partition].Re = 0;
	ptube[partition].Nu = 0;
	ptube[partition].hi = 0;
	ptube[partition].fi = 0;
	ptube[partition].DP = 0;
	
	for(i=0; i<partition; i++) {
	ptube[partition].conductivity+= pcommon[i].length*ptube[i].conductivity;	
	ptube[partition].viscosity += pcommon[i].length*ptube[i].viscosity;
	ptube[partition].density += pcommon[i].length*ptube[i].density;
	ptube[partition].Cp += pcommon[i].length*ptube[i].Cp;
	ptube[partition].Pr += pcommon[i].length*ptube[i].Pr;
	ptube[partition].e_ratio += pcommon[i].length*ptube[i].e_ratio;
	ptube[partition].p_ratio += pcommon[i].length*ptube[i].p_ratio;
	ptube[partition].angle_ratio += pcommon[i].length*ptube[i].angle_ratio;
	ptube[partition].velocity += pcommon[i].length*ptube[i].velocity;
	ptube[partition].Re += pcommon[i].length*ptube[i].Re;
	ptube[partition].Nu += pcommon[i].length*ptube[i].Nu;
	ptube[partition].hi += pcommon[i].length*ptube[i].hi;
	ptube[partition].fi += pcommon[i].length*ptube[i].fi;
	ptube[partition].DP += ptube[i].DP;
	}
	
	ptube[partition].conductivity/= pcommon[partition].length;	
	ptube[partition].viscosity /= pcommon[partition].length;
	ptube[partition].density /= pcommon[partition].length;
	ptube[partition].Cp /= pcommon[partition].length;
	ptube[partition].Pr /= pcommon[partition].length;
	ptube[partition].e_ratio /= pcommon[partition].length;
	ptube[partition].p_ratio /= pcommon[partition].length;
	ptube[partition].angle_ratio /= pcommon[partition].length;
	ptube[partition].velocity /= pcommon[partition].length;
	ptube[partition].Re /= pcommon[partition].length;
	ptube[partition].Nu /= pcommon[partition].length;
	ptube[partition].hi /= pcommon[partition].length;
	ptube[partition].fi /= pcommon[partition].length;
}

void SumNAvg_Annular(Annular* pannular, Common* pcommon) {
	int i;
	pannular[partition].in_temp = 0;	// 무의미한 값들.
	pannular[partition].out_temp = 0;
	pannular[partition].avg_temp = 0;
	pannular[partition].index = 0;
	pannular[partition].in_enthalpy = 0;
	pannular[partition].out_enthalpy = 0;
	
	pannular[partition].conductivity= 0;	//conductivity 부터 시작.
	pannular[partition].viscosity = 0;
	pannular[partition].density = 0;
	pannular[partition].Cp = 0;
	pannular[partition].Pr = 0;
	pannular[partition].e_ratio = 0;
	pannular[partition].p_ratio = 0;
	pannular[partition].r_ratio = 0;//
	pannular[partition].Dhyd = 0;//
	pannular[partition].Aeff = 0;//
	pannular[partition].angle_ratio = 0;
	pannular[partition].velocity = 0;
	pannular[partition].Re = 0;
	pannular[partition].Nu = 0;
	pannular[partition].ho = 0;
	pannular[partition].ef = 0;//
	pannular[partition].fo = 0;
	pannular[partition].DP = 0;
	
	for(i=0; i<partition; i++) {
	pannular[partition].conductivity+= pcommon[i].length*pannular[i].conductivity;	
	pannular[partition].viscosity += pcommon[i].length*pannular[i].viscosity;
	pannular[partition].density += pcommon[i].length*pannular[i].density;
	pannular[partition].Cp += pcommon[i].length*pannular[i].Cp;
	pannular[partition].Pr += pcommon[i].length*pannular[i].Pr;
	pannular[partition].e_ratio += pcommon[i].length*pannular[i].e_ratio;
	pannular[partition].p_ratio += pcommon[i].length*pannular[i].p_ratio;
	pannular[partition].r_ratio += pcommon[i].length*pannular[i].r_ratio;//
	pannular[partition].Dhyd += pcommon[i].length*pannular[i].Dhyd;//
	pannular[partition].Aeff += pcommon[i].length*pannular[i].Aeff;//
	pannular[partition].angle_ratio += pcommon[i].length*pannular[i].angle_ratio;
	pannular[partition].velocity += pcommon[i].length*pannular[i].velocity;
	pannular[partition].Re += pcommon[i].length*pannular[i].Re;
	pannular[partition].Nu += pcommon[i].length*pannular[i].Nu;
	pannular[partition].ho += pcommon[i].length*pannular[i].ho;
	pannular[partition].ef += pcommon[i].length*pannular[i].ef;//
	pannular[partition].fo += pcommon[i].length*pannular[i].fo;
	pannular[partition].DP += pannular[i].DP;
	}
	
	pannular[partition].conductivity/= pcommon[partition].length;	
	pannular[partition].viscosity /= pcommon[partition].length;
	pannular[partition].density /= pcommon[partition].length;
	pannular[partition].Cp /= pcommon[partition].length;
	pannular[partition].Pr /= pcommon[partition].length;
	pannular[partition].e_ratio /= pcommon[partition].length;
	pannular[partition].p_ratio /= pcommon[partition].length;
	pannular[partition].angle_ratio /= pcommon[partition].length;
	pannular[partition].r_ratio /= pcommon[partition].length;//
	pannular[partition].Dhyd /= pcommon[partition].length;//
	pannular[partition].Aeff /= pcommon[partition].length;//
	pannular[partition].velocity /= pcommon[partition].length;
	pannular[partition].Re /= pcommon[partition].length;
	pannular[partition].Nu /= pcommon[partition].length;
	pannular[partition].ho /= pcommon[partition].length;
	pannular[partition].ef /= pcommon[partition].length;//
	pannular[partition].fo /= pcommon[partition].length;
}
	//<\Evaluate Sum of Avg Values>
//<\Fucntion contents>

