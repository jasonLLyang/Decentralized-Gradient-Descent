import java.util.*;
import java.io.*;
public class DGDSimulation {
    private static final String IOFOLDER="..\\..\\DGD_IO\\";
    private static double val(double[] p, double x) {
        return p[0]+p[1]*x+p[2]*x*x;
    }
    private static double grad(double[] p, double x) {
        return p[1]+p[2]*2*x;
    }
    private static class Test {
        int N; String P; long seed;
        int[][] G;
        double[][] Ps;
        double[] locs0;
        Test(int N, String P, long seed, int[][] G, double[][] Ps, double[] locs0) {
            this.N=N;
            this.P=P;
            this.seed=seed;
            this.G=G;
            this.Ps=Ps;
            this.locs0=locs0;
        }
        public String toString() {
            StringBuilder out=new StringBuilder("N,P,seed="+N+","+P+","+seed+"\n");
            for (int i=0; i<N; i++)
                out.append(Arrays.toString(G[i]))
                        .append(" ").append(Arrays.toString(Ps[i]))
                        .append(" ").append(locs0[i]).append("\n");
            return out.toString();
        }
    }
    private static List<Test> graphs(String fin) throws IOException {
        BufferedReader in=new BufferedReader(new FileReader(fin));
        List<Test> out=new ArrayList<>();
        while (true) {
            String param=in.readLine();
            if (param==null) break;
            int N; String P; long seed; {
                String[] tmp=param.split("=")[1].split(",");
                N=Integer.parseInt(tmp[0]);
                P=tmp[1];
                seed=Long.parseLong(tmp[2]);
            }
            int[][] G;
            G=new int[N][];
            for (int i=0; i<N; i++) {
                String[] l=in.readLine().split(":");
                int v=Integer.parseInt(l[0]);
                String[] adjs=l[1].split(",");
                G[v]=new int[adjs.length];
                for (int j=0; j<adjs.length; j++)
                    G[v][j]=Integer.parseInt(adjs[j]);
            }
            double[][] Ps=new double[N][];
            for (int i=0; i<N; i++) {
                String[] l=in.readLine().split(":");
                int v=Integer.parseInt(l[0]);
                String[] info=l[1].split(",");
                Ps[v]=new double[3];
                for (int j=0; j<3; j++)
                    Ps[v][j]=Double.parseDouble(info[j]);
            }
            double[] locs0=new double[N];
            String[] info=in.readLine().split(",");
            for (int j=0; j<N; j++)
                locs0[j]=Double.parseDouble(info[j]);
            out.add(new Test(N,P,seed,G,Ps,locs0));
        }
        return out;
    }
    private static double trimmedMean(double[] A, int amt) {
        if (A.length<=2*amt) throw new RuntimeException("Not enough numbers for trimmed mean.");
        double[] vals=new double[A.length];
        System.arraycopy(A,0,vals,0,A.length);
        Arrays.sort(vals);
        double out=0;
        for (int i=amt; i<A.length-amt; i++)
            out+=vals[i];
        return out/(A.length-2*amt);
    }
    private static double meanScr(double[] tP, double[] locs) {
        double[] scrs=new double[locs.length];
        for (int i=0; i<locs.length; i++)
            scrs[i]=val(tP,locs[i]);
        return trimmedMean(scrs,0);
    }
    public static void main(String[] args) throws Exception {
        long st=System.currentTimeMillis();
        int A=1, ITER=10000;
        boolean equivocate=false;
        List<Test> data=graphs(IOFOLDER+"input.txt");
        //System.out.println(data);
        PrintWriter out=new PrintWriter(new FileWriter(IOFOLDER+"stats_out_A="+A+(equivocate?"_equivocate":"")+".txt"));
        out.printf("A,ITER=%d,%d%n",A,ITER);
        System.out.printf("A,ITER=%d,%d%n",A,ITER);
        for (Test $:data) {
            {
                String intro="N,P,seed="+$.N+","+$.P+","+$.seed;
                if ($.seed%10==0) System.out.println(intro);
                out.println(intro);
            }
            int N=$.N;
            int[][] G=new int[N+A][]; {
                System.arraycopy($.G,0,G,0,N);
                int[] tmp=new int[N];
                for (int i=0; i<N; i++) tmp[i]=i;
                for (int i=N; i<N+A; i++)
                    G[i]=tmp;
            }
            double[][] Ps=$.Ps;
            double[] locs0=$.locs0;

            double[] tP=new double[3];
            for (int i=0; i<N; i++)
                for (int j=0; j<3; j++)
                    tP[j]+=Ps[i][j];
            for (int j=0; j<3; j++) tP[j]/=N;
            double optloc=-tP[1]/(2*tP[2]), optscr=val(tP,optloc);
            out.printf("optloc,optscr=%f,%f%n",optloc,optscr);

            for (double T:new double[] {0.01,0.005,0.002,0.001}) {
                out.println("T="+T);

                double[] locs=new double[N];
                System.arraycopy(locs0,0,locs,0,N);
                double[] meanloc=new double[ITER+1], meandiff=new double[ITER+1];
                meanloc[0]=trimmedMean(locs,0); meandiff[0]=meanScr(tP,locs);
                for (int round=1; round<=ITER; round++) {
                    int[] amts=new int[N];
                    for (int i=0; i<N+A; i++)
                        for (int j:G[i])
                            amts[j]++;
                    double[][] rvals=new double[N][];
                    for (int i=0; i<N; i++) rvals[i]=new double[amts[i]];
                    int[] idxs=new int[N];
                    for (int i=0; i<N+A; i++)
                        for (int j:G[i])
                            rvals[j][idxs[j]++]=i<N?locs[i]://1000_000.0;
                                                            (j%2==0?1000_000.0:-1000_000.0);
                    double[] nlocs=new double[N];
                    for (int i=0; i<N; i++)
                        nlocs[i]=trimmedMean(rvals[i],A)-T*grad(Ps[i],locs[i]);
                    locs=nlocs;
                    meanloc[round]=trimmedMean(locs,0); meandiff[round]=meanScr(tP,locs)-optscr;
                }

                out.println("meanloc_fin="+meanloc[ITER]);
                out.println("meandiff_fin="+meandiff[ITER]);
            }
        }
        out.close();
        System.out.println("time="+(System.currentTimeMillis()-st));
    }
}
