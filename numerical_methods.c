#include <stdio.h>
#include <math.h>
#define MAX 100
double power(double x, int pow_num);
int fakt(int N);
double f(double *result, double coff[], int pow_arr[], int term_num, double x);
double f_diff(double *result, double coff[], int pow_arr[], int term_num, double x);
void get_poly(int *term_num, double coff[], int pow_arr[]);
void print_poly(int term_num, double coff[],int pow_arr[]);
void get_start(double *start);
void get_final(double *final);
void get_tolerance(double *tolerance);
void get_spuare_matrix(int *N, double matrix[][MAX]);
void print_equation_mtr(int N, double temp_matrix[][MAX], double result_arr[]);
void print_sqr_mtr(int N, double matrix[][MAX]);
void copy_mtr(int N, double matrix[][MAX], double temp_matrix[][MAX]);
void swap(int *a, int *b);
void permutation(int start, int final, int matrix[][MAX], int arr[], int *count);
void gauss_upper_triangle(int N, double matrix[][MAX], double temp_matrix[][MAX], double result_arr[]);
void get_equation(int *N, double matrix[][MAX], double result_arr[]);
void identity_marix(int N, double identity_matrix[MAX][MAX]);
void invert_matrix(int N, double matrix[][MAX], double temp_matrix[][MAX], double identity_matrix[][MAX]);
void upper_triangle(int N, double matrix[][MAX], double temp_matrix[][MAX]);
double determinant(int N, double matrix[][MAX], double temp_matrix[][MAX]);
double forvard_diff(int i, int which_diff, double result_arr[]);
int main(){
	int choice,dd=1, N;
	double coff[MAX], result_arr[MAX];
	double matrix[MAX][MAX];
	double temp_matrix[MAX][MAX];
	double temp_arr[MAX];
	double variable_temp[MAX];
	double result_temp[MAX];
	double identity_matrix[MAX][MAX];
	int perm[MAX][MAX];
	int arr[MAX];
	int pow_arr[MAX];
	int term_num, i, j, k, l, a, n, yes, choice2, iteration, max_i, counter;
	double start,final,tolerance,temp,istart,jfin,max;
	double result1=0, result2=0, x, x_n, x_n2, real_x, x0, h, result=0;
	while (dd==1){
		printf("[0] Quit\n[1] Bisection \n[2] Regula-Falsi \n[3] Newton-Rapshon \n[4] Inverse NxN matrix\n");
		printf("[5] Gauss Elemination\n[6] Gauss Seidal\n[7] Numerical Differentiation\n[8] Simpson's Rule \n[9] Trapezoidaloidal Rule\n[10] Gregory-Newton Interpolation\n");
		printf("Choice: ");	
		scanf("%d",&choice);
		if(choice==0){
			printf("\nProgram terminated.");
			dd=0;
		}
		else if(choice==1){//Bisection					
			get_poly(&term_num, coff, pow_arr);
			print_poly(term_num, coff, pow_arr);
			printf("\n\nrange:");			
			get_start(&start);
			get_final(&final);
			get_tolerance(&tolerance);
			result1=f(&result1, coff, pow_arr, term_num, start);				
			result2=f(&result2, coff, pow_arr, term_num, final);
			if(result1*result2==0){
				if(result1==0){
					printf("root: %lf", start);
				}
				else{
					printf("root: %lf", final);
				}
			}
			else if(result1*result2>0){
				printf("Cannot calculate the root with this method.");
			}
			else{
				printf("\nif you know real root, enter '1': ");
				scanf(" %d",&yes);
				if(yes==1){
					printf("enter real root value: ");
					scanf(" %lf",&real_x);
					printf("real root: %lf\n",real_x);
					i=0;
					x=tolerance;
					while(fabs(real_x-x)>tolerance){
						x=(start+final)/2;//Bisection formula
						result=f(&result, coff, pow_arr, term_num, x);		
						if(result1*result<0){
							result2=result;
							final=x;
						}
						else{
							result1=result;
							start=x;
						}
						i++;
						result=0;						
						printf("\n%d.iteration: %lf",i,x);
						printf("\tnew range: [%lf,%lf]",start,final);
						printf("\ttolerance: %lf",fabs(real_x-x));						
					}
					printf("\n\nRoot: %lf\n",x);
				}
				else{					
					i=0;									
					while((final-start)/power(2,i)>tolerance){
						x=(start+final)/2;//Bisection formula
						result=f(&result, coff, pow_arr, term_num, x);					
						if(result1*result<0){			
							result2=result;
							final=x;										
						}
						else{
							result1=result;
							start=x;
						}																		
						result=0;
						i++;
						printf("\n%d.iteration: %lf",i,x);						
						printf("\tnew range: [%lf,%lf]",start,final);
						printf("\ttolerance: %lf",(final-start)/power(2,i));												
					}
					printf("\n\nRoot: %lf\n",x);									
				}
			}
			printf("\n");
		}
		else if(choice==2){//regula-falsi
			get_poly(&term_num, coff, pow_arr);
			print_poly(term_num, coff, pow_arr);
			printf("\n\nrange:");			
			get_start(&start);
			get_final(&final);
			get_tolerance(&tolerance);
			result1=f(&result1, coff, pow_arr, term_num, start);//calculate start value				
			result2=f(&result2, coff, pow_arr, term_num, final);//calculate final value			
			if(result1*result2==0){
				if(result1==0){
					printf("root: %lf", start);
				}
				else{
					printf("root: %lf", final);
				}
			}
			else if(result1*result2>0){
				printf("Cannot calculate the root with this method.");
			}
			else{
				i=0;																
				while((final-start)/power(2,i)>tolerance){
					x=(start*result2-final*result1)/(result2-result1);//regula-falsi formula
					result=f(&result, coff, pow_arr, term_num, x);					
					if(result1*result<0){			
						result2=result;
						final=x;										
					}
					else{
						result1=result;
						start=x;
					}																						
					i++;
					printf("\n%d.iteration: %lf",i,x);						
					printf("\tnew range: [%lf,%lf]",start,final);										
					printf("\ttolerance: %lf",(final-start)/power(2,i));
					result=0;						
				}
				printf("\n\nRoot: %lf\n",x);
			}
			printf("\n");
		}
		else if(choice==3){//newton-raphson
			get_poly(&term_num, coff, pow_arr);
			print_poly(term_num, coff, pow_arr);			
			get_start(&start);
			get_tolerance(&tolerance);
			i=0;
			x_n=start+2*tolerance;
			x_n2=start;
			while(fabs(start-x_n)>tolerance){	
				start=x_n2;		
				x_n2=start-f(&result, coff, pow_arr, term_num, start)/f_diff(&result1, coff, pow_arr, term_num, start);//newton-raphson formula								
				i++;
				printf("\n%d.iteration: %lf",i,x_n2);
				printf("\ttolerance: %lf",fabs(start-x_n2));
				x_n=x_n2;				
			}
			printf("\n\nRoot: %lf\n",x_n2);
			printf("\n");
		}
		else if(choice==4){//inverse matrix						
			get_spuare_matrix(&N, matrix);
			printf("\nMatrix you entered:\n");
			print_sqr_mtr(N, matrix);
			printf("\n");
			if(determinant(N, matrix, temp_matrix)==0){
				printf("Matrix has no inverse\n");
			}
			else{
				identity_marix(N, identity_matrix);
				invert_matrix(N, matrix, temp_matrix, identity_matrix);
				printf("Inverse matrixi:\n");								
				print_sqr_mtr(N, identity_matrix);
			}			
			printf("\n");
		}
		else if(choice==5){//Gauss Elemination			
			get_equation(&N, matrix, result_arr);
			printf("\nMatrix you entered:\n");
			print_equation_mtr(N, matrix, result_arr);
			printf("\n");
			if(determinant(N,matrix,temp_matrix)==0){
				printf("\nCannot solve it because determinant of the matrix you entered is 0.\n");
			}
			else{
				gauss_upper_triangle(N, matrix, temp_matrix, result_arr);
				for(i=0; i<N; i++){
					temp=temp_matrix[i][i];
					for(j=0; j<N; j++){
						temp_matrix[i][j]=temp_matrix[i][j]/temp;
					}
					result_arr[i]=result_arr[i]/temp;
				}											
				for(i=(N-2); i>=0; i--){
					a=N-1;
					for(j=i; j<=(N-2); j++){
						result_arr[i]=result_arr[i]-result_arr[a]*temp_matrix[i][a];
						a--;
					}						
				}
				printf("\n");			
				for(i=0; i<N; i++){
					printf("value of %d.variable: %lf\n",i,result_arr[i]);
				}
			}
						
			printf("\n");
		}
		else if(choice==6){//Gauss Seidal
			get_equation(&N, matrix, result_arr);
			printf("\n");
			print_equation_mtr(N, matrix, result_arr);
			printf("\n");						
			counter=0;
			for(i=0; i<N; i++){
				arr[i]=i;
			}
			permutation(0,N-1,perm,arr, &counter);
			max=0;
			for(k=0; k<counter; k++){
				result=1;
				for(i=0;i<N;i++){
					for(j=0;j<N;j++){
						temp_matrix[i][j]=matrix[perm[k][i]][j];
					}
				}
				for(i=0;i<N;i++){
					result=result*fabs(temp_matrix[i][i]);
				}
				if(result>max){
					max=result;
					max_i=k;
				}
			}
			for(i=0;i<N;i++){
				for(j=0;j<N;j++){
					temp_matrix[i][j]=matrix[perm[max_i][i]][j];
				}
				result_temp[i]=result_arr[perm[max_i][i]];
			}											
			printf("diagonally dominant matrix:\n");
			print_equation_mtr(N, temp_matrix, result_temp);
			printf("\n");
			for(i=0; i<N; i++){
				printf("Enter the start value of %d.variable:",i);
				scanf(" %lf", &temp_arr[i]);
			}
			
			printf("\nStopping conditions:\n[1]iteration value\n[2]tolerance value\nchoice: ");
			scanf(" %d",&choice2);
			if(choice2==1){
				printf("\nEnter a max iteration value: ");
				scanf(" %d", &iteration);			
				printf("\nstart values of variables:\n");
				for(i=0; i<N; i++){
					printf("%d.variable: %lf\t",i,temp_arr[i]);
				}
				printf("\n\n");
				for(k=0; k<iteration; k++){	
					for(i=0; i<N; i++){
						variable_temp[i]=temp_arr[i];
					}								
					for(i=0; i<N; i++){
						result=result_temp[i];
						for(j=0; j<N; j++){
							if(j!=i){						
								result=result-temp_matrix[i][j]*temp_arr[j];						
							}
						}
						temp_arr[i]=result/temp_matrix[i][i];					
					}
					printf("%d. iteration:\n",k);
					for(i=0; i<N; i++){
						printf("%d.variable: %lf, \t",i,temp_arr[i]);
						variable_temp[i]=fabs(temp_arr[i]-variable_temp[i]);
						printf("tolerance: %lf\n",variable_temp[i]);
					}
					printf("\n\n");
				}
				printf("Roots:\n");
				for(i=0;i<N;i++){
					printf("%d.Root: %lf\n",i,temp_arr[i]);
				}	
				printf("\n");
			}
			else if(choice2==2){
				get_tolerance(&tolerance);
				printf("\nstart values of variables:\n");
				for(i=0; i<N; i++){
					printf("%d.variable: %lf\t",i,temp_arr[i]);
				}
				printf("\n\n");
				counter=0;
				k=0;
				while(counter<N){
					counter=0;								
					for(i=0; i<N; i++){
						variable_temp[i]=temp_arr[i];
					}													
					for(i=0; i<N; i++){
						result=result_temp[i];
						for(j=0; j<N; j++){
							if(j!=i){						
								result=result-temp_matrix[i][j]*temp_arr[j];						
							}
						}
						temp_arr[i]=result/temp_matrix[i][i];					
					}
					printf("%d. iteration:\n",k);
					for(i=0; i<N; i++){
						printf("%d.variable: %lf, \t",i,temp_arr[i]);
						variable_temp[i]=fabs(temp_arr[i]-variable_temp[i]);
						printf("tolerance: %lf\n",variable_temp[i]);
					}
					k++;
					printf("\n\n");
					for(i=0; i<N; i++){
						if(variable_temp[i]<tolerance){
							counter++;
						}					
					}	
				}
				printf("Roots:\n");
				for(i=0;i<N;i++){
					printf("%d.Root: %lf\n",i,temp_arr[i]);
				}									
			}																	
			printf("\n");
		}
		else if(choice==7){//Numerical Differentiation
			printf("\n[1]Backward Difference\n[2]Forward Difference\n[3]Centered Difference\nChoice: ");
			scanf(" %d", &choice2);
			get_poly(&term_num, coff, pow_arr);
			print_poly(term_num, coff, pow_arr);
			printf("\nEnter the x value: ");
			scanf(" %lf", &x);
			printf("Enter the h value: ");
			scanf(" %lf", &h);
			if(h==0){
				printf("Cannot the calculate it because the h value is 0.");
			}
			else{
				if(choice2==1){
					printf("result: %lf", (f(&result, coff, pow_arr, term_num, x)-f(&result, coff, pow_arr, term_num, (x-h)))/h);
				}
				else if(choice2==2){
					printf("result: %lf", (f(&result, coff, pow_arr, term_num, (x+h))-f(&result, coff, pow_arr, term_num, x))/h);
				}
				else if(choice2==3){
					printf("result: %lf", (f(&result, coff, pow_arr, term_num, (x+h))-f(&result, coff, pow_arr, term_num, (x-h)))/(2*h));
				}
				else{
					printf("Wrong choice!");
				}
			}			
			printf("\n");
		}
		else if(choice==8){//Simpson's Rule
			printf("\n[1]Simpson 1/3\n[2]Simpson 3/8\nChoice: ");
			scanf(" %d", &choice2);
			get_poly(&term_num, coff, pow_arr);
			print_poly(term_num, coff, pow_arr);
			printf("\n\nintegral's");
			get_start(&start);
			get_final(&final);			
			printf("Enter the n value: ");
			scanf(" %d",&n);			
			h=(final-start)/n;
			result=0;		
			if(choice2==1){//simpson 1/3
				for(istart=start; istart<final; istart+=h){//simpson 1/3 formula
					jfin=istart+h;					
					result+=(jfin-istart)*(f(&result1, coff, pow_arr, term_num, istart)+4*f(&result1, coff, pow_arr, term_num, ((istart+jfin)/2))+f(&result1, coff, pow_arr, term_num, jfin))/6;
				}
				printf("\nresult: %lf",result);
			}
			else if(choice2==2){//simpson 3/8
				for(istart=start; istart<final; istart+=h){//simpson 3/8 formula
					jfin=istart+h;
					result+=(jfin-istart)*(f(&result1, coff, pow_arr, term_num, istart)+3*f(&result1, coff, pow_arr, term_num, istart+(jfin-istart)/3)+3*f(&result1, coff, pow_arr, term_num, istart+2*(jfin-istart)/3)+f(&result1, coff, pow_arr, term_num, jfin))/8;
				}
				printf("\nresult: %lf",result);
			}
			else{
				printf("Wrong choice!");
			}			
			printf("\n");
		}
		else if(choice==9){//Trapezoidal 
			get_poly(&term_num, coff, pow_arr);
			print_poly(term_num, coff, pow_arr);
			printf("\n\nintegral's");
			get_start(&start);
			get_final(&final);
			printf("Enter the n value: ");//
			scanf(" %d",&n);
			h=(final-start)/n;
			for(istart=start; istart<final; istart+=h){//Trapezoidal formula
				jfin=istart+h;
				result+=(jfin-istart)*(f(&result1, coff, pow_arr, term_num, istart)+f(&result1, coff, pow_arr, term_num, jfin))/2;
			}
			printf("\nresult: %lf",result);
			printf("\n");
		}
		else if(choice==10){//Gregory Newton Enterpolation
			printf("\nEnter the start value of x (x0): ");
			scanf(" %lf", &x0);
			printf("Enter the difference between of x's (h value): ");
			scanf(" %lf", &h);
			printf("Enter the counts of f(x): ");
			scanf(" %d", &n);
			printf("Enter the x value that you want to calculate: ");
			scanf(" %lf", &x);
			temp=x0;
			for(i=0; i<n; i++){
				printf("Enter the f(%lf) value: ", temp);
				scanf(" %lf", &result_arr[i]);
				temp=temp+h;
			}
			result2=0;
			for(i=1; i<n; i++){
				temp=x0;
				result=1;
				for(j=1; j<=i; j++){
					result=result*(x-temp)/(h*j);
					temp=temp+h;
				}
				result=result*forvard_diff(0,i,result_arr);
				result2=result2+result;
			}
			result2 +=result_arr[0];
			printf("\nresult: %lf", result2);
			
			printf("\n");
		}
		else{
			printf("\nWrong access!\n\n");
		}
	}
	return 0;
}
double power(double x, int pow_num){
	int i;
	double result=1;
	for(i=0;i<pow_num;i++){
		result=result*x;
	}
	return result;
}
int fakt(int N){//factorial
	int result=1, i;
	for(i=1; i<=N; i++){
		result=result*i;
	}
	return result;
}
double f(double *result, double coff[], int pow_arr[], int term_num, double x){//calculate the polynomial function
	int i;
	*result=0;
	for(i=0;i<term_num;i++){
		*result=*result+coff[i]*power(x,pow_arr[i]);
	}
	return *result;
}
double f_diff(double *result, double coff[], int pow_arr[], int term_num, double x){//calculate the differentiation of polynomial function
	int i;
	*result=0;
	for(i=0;i<term_num;i++){
		if(pow_arr[i]!=0){
			*result=*result+coff[i]*pow_arr[i]*power(x,(pow_arr[i]-1));
		}				
	}
	if(result!=0){
		return *result;
	}
	else{
		return 0.000001;
	}	
}
void get_poly(int *term_num, double coff[], int pow_arr[]){
	int i,dd;
	printf("\nEnter the count of term of polynomial: ");
	scanf("%d",term_num);
	for(i=0;i<*term_num;i++){		
		printf("Enter the coefficient of the %d.term: ",i);
		scanf("%lf",&coff[i]);
		dd=1;
		while(dd==1){
			printf("Enter the power of the %d.term: ",i);
			scanf("%d",&pow_arr[i]);
			if(pow_arr[i]<0){
				printf("Power cannot be negative!\n");
				dd=1;
			}	
			else{
				dd=0;
			}
		}				
	}	
}
void print_poly(int term_num, double coff[],int pow_arr[]){
	int i;
	printf("\nPolynomial: ");
	for(i=0;i<(term_num-1);i++){
		printf("%lfx^%d + ",coff[i],pow_arr[i]);
	}
	printf("%lfx^%d ",coff[i],pow_arr[i]);
}
void get_start(double *start){
	printf("\nstart value: ");
	scanf(" %lf", start);	
}
void get_final(double *final){
	printf("final value: ");
	scanf(" %lf", final);
}	
void get_tolerance(double *tolerance){
	printf("Enter the tolerance value: ");
	scanf(" %lf",tolerance);
}
void get_spuare_matrix(int *N, double matrix[][MAX]){
	int i,j;
	printf("Enter the size of square matrix: ");
	scanf(" %d",N);
	for(i=0;i<*N;i++){
		for(j=0;j<*N;j++){
			printf("Enter the %d,%d.element of matrix: ",i,j);
			scanf(" %lf",&matrix[i][j]);
		}
	}
}
void print_equation_mtr(int N, double temp_matrix[][MAX], double result_arr[]){
	int i, j;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%lf ", temp_matrix[i][j]);
		}
		printf("| %lf",result_arr[i]);
		printf("\n");
	}
}
void print_sqr_mtr(int N, double matrix[][MAX]){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("%lf\t",matrix[i][j]);
		}
		printf("\n");
	}
}
void copy_mtr(int N, double matrix[][MAX], double temp_matrix[][MAX]){
	int i,j;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			temp_matrix[i][j]=matrix[i][j];
		}
	}
}
void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}
void permutation(int start, int final, int matrix[][MAX], int arr[], int *count){
	int i;
	if(start==final){
		for(i=0;i<=final;i++){
			matrix[*count][i]=arr[i];
		}
		(*count)++;
		return;
	}
	for(i=start; i<=final; i++){
		swap(&arr[start],&arr[i]);
		permutation(start+1,final, matrix, arr, count);
		swap(&arr[start],&arr[i]);
	}
}
void gauss_upper_triangle(int N, double matrix[][MAX], double temp_matrix[][MAX], double result_arr[]){
	int i,j,k,l,a;
	double temp;
	copy_mtr(N, matrix, temp_matrix);
	a=1;
	for(k=0;k<(N-1); k++){	
		l=k;		
		for(i=a;i<N;i++){
			while(l<(N-1) && temp_matrix[l][k]==0){
				l++;				
			}
			while(k<(N-1) && temp_matrix[l][k]==0){
				k++;
			}
			temp=temp_matrix[i][k]/temp_matrix[l][k];			
			for(j=0;j<N;j++){				
				temp_matrix[i][j]=temp_matrix[i][j]-temp_matrix[l][j]*temp;				
			}
			result_arr[i]=result_arr[i]-result_arr[l]*temp;
		}
		a++;				
	}	
}
void get_equation(int *N, double matrix[][MAX], double result_arr[]){
	int i, j;
	printf("\nEnter the count of equation or variable: ");
	scanf(" %d",N);
	for(i=0; i<*N; i++){
		printf("%d.equation's;\n ",i);
		for(j=0; j<*N; j++){
			printf("%d.variable' coefficient: ",j);
			scanf(" %lf", &matrix[i][j]);
		}
		printf("%d.equation's result: ",i);
		scanf(" %lf", &result_arr[i]);
	}
}
void identity_marix(int N, double identity_matrix[][MAX]){
	int i,j;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			if(i==j){
				identity_matrix[i][j]=1;
			}
			else{
				identity_matrix[i][j]=0;
			}			
		}
	}
}
void invert_matrix(int N, double matrix[][MAX], double temp_matrix[][MAX], double identity_matrix[][MAX]){
	int i,j,k,l,a;
	double temp, temp2;
	copy_mtr(N, matrix, temp_matrix);
	a=1;
	for(k=0;k<(N-1); k++){	
		l=k;		
		for(i=a;i<N;i++){
			while(l<(N-1) && temp_matrix[l][k]==0){
				l++;				
			}
			while(k<(N-1) && temp_matrix[l][k]==0){
				k++;
			}
			temp=temp_matrix[i][k]/temp_matrix[l][k];			
			for(j=0;j<N;j++){				
				temp_matrix[i][j]=temp_matrix[i][j]-temp_matrix[l][j]*temp;	
				identity_matrix[i][j]=identity_matrix[i][j]-identity_matrix[l][j]*temp;			
			}			
		}
		a++;				
	}
	for(i=0; i<N; i++){
		temp=temp_matrix[i][i];
		for(j=0; j<N; j++){
			temp_matrix[i][j]=temp_matrix[i][j]/temp;
			identity_matrix[i][j]=identity_matrix[i][j]/temp;
		}		
	}		
	for(i=(N-2); i>=0; i--){
		a=N-1;	
		for(l=(N-1); l>i; l--){
			temp=temp_matrix[i][l];	
			for(j=(N-1); j>=0; j--){										
				temp_matrix[i][j]=temp_matrix[i][j]-temp_matrix[a][j]*temp;
				identity_matrix[i][j]=identity_matrix[i][j]-identity_matrix[a][j]*temp;					
			}
			a--;	
		}				
	}		
}
void upper_triangle(int N, double matrix[][MAX], double temp_matrix[][MAX]){
	int i,j,k,l,a;
	double temp;
	copy_mtr(N, matrix, temp_matrix);
	a=1;
	for(k=0;k<(N-1); k++){	
		l=k;		
		for(i=a;i<N;i++){
			while(l<(N-1) && temp_matrix[l][k]==0){
				l++;				
			}
			while(k<(N-1) && temp_matrix[l][k]==0){
				k++;
			}
			temp=temp_matrix[i][k]/temp_matrix[l][k];			
			for(j=0;j<N;j++){				
				temp_matrix[i][j]=temp_matrix[i][j]-temp_matrix[l][j]*temp;				
			}			
		}
		a++;				
	}
}
double determinant(int N, double matrix[][MAX], double temp_matrix[][MAX]){
	int det=1, i;
	upper_triangle(N, matrix, temp_matrix);
	for(i=0; i<N; i++){
		det=det*temp_matrix[i][i];
	}
	return det;	
}
double forvard_diff(int i, int which_diff, double result_arr[]){
	if(which_diff>0){
		return forvard_diff(i+1,which_diff-1,result_arr)-forvard_diff(i,which_diff-1,result_arr);
	}
	else{
		return result_arr[i];
	}	
}
