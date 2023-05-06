import java.util.Random;
import java.awt.Graphics;
import java.awt.Color;
import java.util.Arrays;
import java.util.ArrayList;

public final class lenia {
    int height;
    int width;
    double[][] state;
    double[][] mem;
    int R;
    int T;
    double dx;
    double dt;
    double mu;
    double sigma;
    double alpha;
    double[] beta;
    int B;
    double[][] ker;
    
    int plot_size;
    Color[] arr_clr;
    Random rand=new Random();
    int[][][] coord;
    ArrayList<int[][]> coordList;

    lenia(double[] beta,double mu,double sigma,int R,int T,int h,int w,int ps){
        this.T=T;
        this.dt=1.0/T;
        this.R=R;
        this.dx=1.0/R;
        this.alpha=4.0;
        this.beta=beta;
        this.B=beta.length;
        height = h;
        width = w;
        plot_size=ps;
        this.mu=mu;this.sigma=sigma;
        state=zero_arr_d(height,width);
        mem=zero_arr_d(height,width);
        ker=make_ker_exp(this.alpha,this.beta,this.B);
        arr_clr=new Color[256];
        for(int i=0;i<256;i++){
            arr_clr[i]=new Color(i,i,i);
        }
        coord=make_coord(height,width);
        coordList=new ArrayList<int[][]>(Arrays.asList(coord));
        //coordStream=coordList.parallelStream();
    }
    int[][][] make_coord(int h,int w){
        int[][][] cd=new int[h][w][2];
        for (int i=0;i<h;i++){
            for(int j=0;j<w;j++){
                cd[i][j][0]=i;
                cd[i][j][1]=j;
            }
        }
        return cd;
    }

    double[][] zero_arr_d(int h,int w){
        double[][] arr=new double[h][w];
        for(int i=0;i<h;i++){
            for(int j=0;j<w;j++){
                arr[i][j]=0;
            }
        }
        return arr;
    }

    void rand_cir_arr(int num_cir){
        int r;
        int x;
        int y;
        int dx;
        int dy;
        int[] pos;
        for(int i=0;i<num_cir;i++){
            r=rand.nextInt(7,13);
            x=rand.nextInt(0,width);
            y=rand.nextInt(0,height);
            for(int j=0;j<=2*r+1;j++){
                for(int k=0;k<=2*r+1;k++){
                    dx=k-r;dy=j-r;
                    if(dx*dx+dy*dy<=r*r){
                        pos=periodic_pos(x,y,dx,dy);
                        state[pos[1]][pos[0]]=1.0;
                    }
                }
            }
        }
    }

    double[][] make_ker_exp(double alpha,double[] beta,int B){
        double ux;
        double uy;
        double r=0;
        double[][] k=zero_arr_d(2*R+1,2*R+1);
        double s=0;
        for(int i=0;i<2*R+1;i++){
            for(int j=0;j<2*R+1;j++){
                ux=(-R+j)*dx;
                uy=(-R+i)*dx;
                r=Math.pow(ux*ux+uy*uy,0.5);
                //System.out.println(r);
                if(r<1.0){
                    double br_a=(B*r)%(1.0);
                    int br_i=(int)(Math.floor(B*r));
                    k[i][j]=beta[br_i]*Math.exp(alpha-alpha/(4*br_a*(1-br_a)+0.000001));
                    s+=k[i][j];
                }
            }
        }
        for(int i=0;i<2*R+1;i++){
            for(int j=0;j<2*R+1;j++){
                k[i][j]=k[i][j]/(s+0.000001);
            }
        }
        return k;
    }

    int[] periodic_pos(int x,int y,int dx,int dy){
        int h=this.height;
        int w=this.width;
        int nx=x+dx;
        int ny=y+dy;
        if(nx<0){nx=nx+w;}
        else if(nx>w-1){nx=nx-w;}
        if(ny<0){ny=ny+h;}
        else if(ny>h-1){ny=ny-h;}
        int[] pos=new int[]{nx,ny};
        return pos;

    }

    double calc_conv_p(double[][] arr,double[][] ker,int rk,int x,int y){
        double cv=0;
        int[] pos;
        for(int i=0;i<rk*2+1;i++){
            for(int j=0;j<rk*2+1;j++){
                pos=periodic_pos(x,y,j-rk,i-rk);
                cv=cv+ker[i][j]*arr[pos[1]][pos[0]];
            }
        }
        return cv;
    }

    void rand_state(){
        int t=5;
        for (int i=0;i<height;i++){
            for(int j=0;j<width;j++){
                if(i>height/t&&i<height/t*2
                &&j>width/t&&j<width/t*2){
                state[i][j]=Math.pow(
                    rand.nextDouble(0.0,1.0),2);
                }
            }
        }
    }

    double growth_bell(double u){
        return 2*Math.exp(-(u-mu)*(u-mu)/(2*sigma*sigma))-1;
    }

    double clip(double s,double mini,double maxx){
        if (s<mini){return mini;}
        else if(s>maxx){return maxx;}
        else{return s;}
    }

    void calc_update(){
        int h=this.height;
        int w=this.width;
        //List<int[][]> clis = new ArrayList<int[][]>(Arrays.asList(coord));
        coordList.parallelStream().forEach(
            arr->{
                double u;
                double g;
                double s;
                for(int j=0;j<w;j++){
                    u=calc_conv_p(state,ker,R,arr[j][1],arr[j][0]);
                    g=growth_bell(u);
                    s=clip(state[arr[j][0]][arr[j][1]]+g*dt, 0.0, 1.0);
                    mem[arr[j][0]][arr[j][1]]=s;//state[i][j]+dt*(2*S-1);
                    //if(S<-1){printd(S);}
                    if(mem[arr[j][0]][arr[j][1]]<0){
                        mem[arr[j][0]][arr[j][1]]=0;
                    }else if(mem[arr[j][0]][arr[j][1]]>1){
                        mem[arr[j][0]][arr[j][1]]=1;
                    } 
                }
            }
        );
        
    }

    void state_update(){
        for(int i=0;i<height;i++){
            for(int j=0;j<width;j++){
                state[i][j]=mem[i][j];
            }
        }
    }

    void printd(double x){
        System.out.printf("%f\n",x);
    }

    void disp_life(Graphics g){
        int m=plot_size;

        double c;
        for(int i=0;i<height;i++){
            for(int j=0;j<width;j++){
                c=255*state[i][j];
                g.setColor(arr_clr[(int)c]);
                g.fillRect(j*m,i*m,m,m);
            }
        }
        //plot_kernel(g, ker_r,ra);
    }

    void plot_kernel(Graphics g,double[][] ker,int r){
        double c;
        int m=plot_size;
        double d=0;
        for(int i=0;i<2*r+1;i++){
            for(int j=0;j<2*r+1;j++){
                c=255*ker[i][j];
                d+=ker[i][j];
                System.out.printf("%f ",ker[i][j]);
                g.setColor(arr_clr[(int)c]);
                g.fillRect(j*m*2,i*m*2,m*2,m*2);
            }
            System.out.printf("\n");
        }
        printd(d);
    }

    void make_orbium_at(int y, int x){
        double[][] orb_arr=new double[][]{
            {0,0,0,0,0,0,0.1,0.14,0.1,0,0,0.03,0.03,0,0,0.3,0,0,0,0},
            {0,0,0,0,0,0.08,0.24,0.3,0.3,0.18,0.14,0.15,0.16,0.15,0.09,0.2,0,0,0,0},
            {0,0,0,0,0,0.15,0.34,0.44,0.46,0.38,0.18,0.14,0.11,0.13,0.19,0.18,0.45,0,0,0},
            {0,0,0,0,0.06,0.13,0.39,0.5,0.5,0.37,0.06,0,0,0,0.02,0.16,0.68,0,0,0},
            {0,0,0,0.11,0.17,0.17,0.33,0.4,0.38,0.28,0.14,0,0,0,0,0,0.18,0.42,0,0},
            {0,0,0.09,0.18,0.13,0.06,0.08,0.26,0.32,0.32,0.27,0,0,0,0,0,0,0.82,0,0},
            {0.27,0,0.16,0.12,0,0,0,0.25,0.38,0.44,0.45,0.34,0,0,0,0,0,0.22,0.17,0},
            {0,0.07,0.2,0.02,0,0,0,0.31,0.48,0.57,0.6,0.57,0,0,0,0,0,0,0.49,0},
            {0,0.59,0.19,0,0,0,0,0.2,0.57,0.69,0.76,0.76,0.49,0,0,0,0,0,0.36,0},
            {0,0.58,0.19,0,0,0,0,0,0.67,0.83,0.9,0.92,0.87,0.12,0,0,0,0,0.22,0.07},
            {0,0,0.46,0,0,0,0,0,0.7,0.93,1,1,1,0.61,0,0,0,0,0.18,0.11},
            {0,0,0.82,0,0,0,0,0,0.47,1,1,0.98,1,0.96,0.27,0,0,0,0.19,0.1},
            {0,0,0.46,0,0,0,0,0,0.25,1,1,0.84,0.92,0.97,0.54,0.14,0.04,0.1,0.21,0.05},
            {0,0,0,0.4,0,0,0,0,0.09,0.8,1,0.82,0.8,0.85,0.63,0.31,0.18,0.19,0.2,0.01},
            {0,0,0,0.36,0.1,0,0,0,0.05,0.54,0.86,0.79,0.74,0.72,0.6,0.39,0.28,0.24,0.13,0},
            {0,0,0,0.01,0.3,0.07,0,0,0.08,0.36,0.64,0.7,0.64,0.6,0.51,0.39,0.29,0.19,0.04,0},
            {0,0,0,0,0.1,0.24,0.14,0.1,0.15,0.29,0.45,0.53,0.52,0.46,0.4,0.31,0.21,0.08,0,0},
            {0,0,0,0,0,0.08,0.21,0.21,0.22,0.29,0.36,0.39,0.37,0.33,0.26,0.18,0.09,0,0,0},
            {0,0,0,0,0,0,0.03,0.13,0.19,0.22,0.24,0.24,0.23,0.18,0.13,0.05,0,0,0,0},
            {0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.09,0.07,0.05,0.01,0,0,0,0,0}
        };
        for(int i=0;i<40;i++){
            for(int j=0;j<40;j++){
                int[] pos=periodic_pos(x,y,j,i);
                state[pos[1]][pos[0]]=orb_arr[i/2][j/2];
            }
        }
    }
}
