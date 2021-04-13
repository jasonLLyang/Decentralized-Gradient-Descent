import java.util.*;
public class DGD {
    //for now only focus on the 1D case
    //all f_i(x) are s.t. f_i'(x) are L-Lipschitz continuous
    private static double L=2, gradstep=0.0001;
    private static double eval(double[] coeffs, double x) {
        double out=0, px=1;
        for (int i=0; i<coeffs.length; i++, px*=x)
            out+=coeffs[i]*px;
        return out;
    }
    private static double grad(double[] coeffs, double x) {
        return coeffs[1]+2*coeffs[2]*x;//(eval(coeffs,x+gradstep)-eval(coeffs,x))/gradstep;
    }
    public static void main(String[] args) {
        SplittableRandom rnd=new SplittableRandom(1);
        int N=100;
        double[][] F=new double[N][];
        for (int i=0; i<N; i++)
            F[i]=new double[] {rnd.nextDouble()*100,rnd.nextDouble()*10,rnd.nextDouble()};
        double[] locs=new double[N];
        for (int i=0; i<N; i++)
            locs[i]=rnd.nextDouble()*50;
        double t=0.5/L;
        /**
         * SUGGESTIONS:
         * random non-fully-connected graph
         * more advanced advesary
         */
        for (int rep=0; rep<200; rep++) {
            double[] tmp=new double[N];
            //for now assume everyone uses exact median
            System.arraycopy(locs,0,tmp,0,N);
            //apply adversary's actions
            for (int i=0; i<40; i++)
                tmp[i]=1000000;
            Arrays.sort(tmp);
            /*double med=0;
            int iamt=N/10;
            for (int i=iamt; i<N-iamt; i++)
                med+=tmp[i];
            med/=N-2*iamt;*/
            double med=tmp[N/2];
            {
                double scr=0;
                for (int i=0; i<N; i++)
                    scr+=eval(F[i],med);
                //if (rep%100==0)
                    System.out.println(scr+" "+med);//+" "+Arrays.toString(locs));
            }
            for (int i=0; i<N; i++)
                locs[i]=med-t*grad(F[i],med);
        }
        {
            double[] tot=new double[3];
            for (int i=0; i<N; i++)
                for (int j=0; j<3; j++)
                    tot[j]+=F[i][j];
            double optx=-tot[1]/(2*tot[2]);
            System.out.println("opt="+eval(tot,optx)+" optx="+optx);//+" "+Arrays.toString(tot));
        }
    }
}
