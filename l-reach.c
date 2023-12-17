#include <stdio.h>
#include <stdlib.h>


int main(){
    int i,j,count=0,c1;
    char c[5]={0,1,2,3,4};
    char d[1000][5]={0};
    

    for(i=0;i<5;i++){
        d[count][0]=c[i];
        d[count][1]=-1;
        d[count++][2]=-1;
        //count++;
        for(j=0;j<5;j++){
            if(i==j){
            d[count][1]=c[j];
            d[count][2]=-2;
            }
            if(i!=j){
            d[count][0]=c[i];
            for(c1=count;c1<count+5;c1++){
            d[c1][1]=c[j];
            //if((c1+1==)
            printf("C%d,",c1);
            }
            //count++;
            printf("b%d,",c[j]);
            for(int k=0;k<5;k++){
                if(j==k){
                d[count][0]=c[i];
                d[count][1]=c[j];
                d[count++][2]=-3;
                }
                //count++;
                if(j!=k){
                printf("c%d,",c[k]);
                d[count][0]=c[i];
                d[count][1]=c[j];
                d[count++][2]=c[k];
                //count++;
                }
            }
            printf("count=%d\n",count);
            }
        }
    }
    
        for(i=0;i<count;i++){
            printf("%d:",i);
        for(j=0;j<3;j++){
        //if(i!=0 && i!=1 && i%5!=1)
        printf("%d,",d[i][j]);
        }
        printf("\n");
    }

    return 0;
}