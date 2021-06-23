import java.util.*;
import java.io.*;
public class DGDSimulation {
    private static final String IOFOLDER="..\\..\\DGD_IO\\"; //<--modify file directory if needed
    private static double val(double[] f, double x) {
        return f[0]+f[1]*x+f[2]*x*x;
    }
    private static double grad(double[] f, double x) {
        return f[1]+f[2]*2*x;
    }
    private static class Test {
        int N; String P; long seed;
        int[][] G;
        double[][] funcs;
        double[] locs0;
        Test(int N, String P, long seed, int[][] G, double[][] funcs, double[] locs0) {
            this.N=N;
            this.P=P;
            this.seed=seed;
            this.G=G;
            this.funcs=funcs;
            this.locs0=locs0;
        }
        public String toString() {
            StringBuilder out=new StringBuilder("N,P,seed="+N+","+P+","+seed+"\n");
            for (int i=0; i<N; i++)
                out.append(Arrays.toString(G[i]))
                        .append(" ").append(Arrays.toString(funcs[i]))
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
            double[][] funcs=new double[N][];
            for (int i=0; i<N; i++) {
                String[] l=in.readLine().split(":");
                int v=Integer.parseInt(l[0]);
                String[] info=l[1].split(",");
                funcs[v]=new double[3];
                for (int j=0; j<3; j++)
                    funcs[v][j]=Double.parseDouble(info[j]);
            }
            double[] locs0=new double[N];
            String[] info=in.readLine().split(",");
            for (int j=0; j<N; j++)
                locs0[j]=Double.parseDouble(info[j]);
            out.add(new Test(N,P,seed,G,funcs,locs0));
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
    private static double meanScr(double[] tfunc, double[] locs) {
        double[] scrs=new double[locs.length];
        for (int i=0; i<locs.length; i++)
            scrs[i]=val(tfunc,locs[i]);
        return trimmedMean(scrs,0);
    }
    public static void main(String[] args) throws Exception {
        long st=System.currentTimeMillis();
        boolean directed=false;
        int A=1, ITER=10000;
        boolean equivocate=true;
        List<Test> data=graphs(IOFOLDER+"input"+(directed?"_directed":"")+".txt");
        System.out.println("# tests="+data.size());
        //System.out.println(data);
        PrintWriter out=new PrintWriter(new FileWriter(IOFOLDER+"out"+(directed?"_directed":"")+"_A="+A+(A>0&&equivocate?"_equivocate":"")+".txt"));
        out.printf("A,ITER=%d,%d%n",A,ITER);
        System.out.printf("A,ITER=%d,%d%n",A,ITER);
        for (Test $:data) {
            {
                String intro="N,P,seed="+$.N+","+$.P+","+$.seed;
                if ($.seed%10==0) System.out.println(intro+" time="+(System.currentTimeMillis()-st));
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
            int[] indeg=new int[N];
            for (int i=0; i<N+A; i++)
                for (int j:G[i])
                    indeg[j]++;
            double[][] funcs=$.funcs;
            double[] locs0=$.locs0;

            double[] tfunc=new double[3];
            for (int i=0; i<N; i++)
                for (int j=0; j<3; j++)
                    tfunc[j]+=funcs[i][j];
            for (int j=0; j<3; j++) tfunc[j]/=N;
            final double optloc=-tfunc[1]/(2*tfunc[2]), optscr=val(tfunc,optloc);
            out.printf("optloc,optscr=%f,%f%n",optloc,optscr);

            for (double T:new double[] {0.01,0.005,0.002,0.001}) {
                out.println("T="+T);
                double[] locs=new double[N]; System.arraycopy(locs0,0,locs,0,N);
                double[] msoldiff=new double[ITER+1], mobjdiff=new double[ITER+1];
                msoldiff[0]=trimmedMean(locs,0)-optloc; mobjdiff[0]=meanScr(tfunc,locs)-optscr;
                for (int round=1; round<=ITER; round++) {
                    double[][] rvals=new double[N][];
                    for (int i=0; i<N; i++) rvals[i]=new double[indeg[i]];
                    int[] idxs=new int[N];
                    for (int i=0; i<N+A; i++)
                        for (int j:G[i])
                            rvals[j][idxs[j]++]=i<N?locs[i]:
                                    (equivocate?(j<N/2?1000_000.0:-1000_000.0):
                                            1000_000.0);
                    /*for (int i=0; i<N; i++)
                        if (idxs[i]!=indeg[i])
                            throw new RuntimeException(idxs[i]+" "+indeg[i]);*/
                    double[] nlocs=new double[N];
                    for (int i=0; i<N; i++)
                        nlocs[i]=(rvals[i].length<=2*A?locs[i]:trimmedMean(rvals[i],A))-T*grad(funcs[i],locs[i]);
                    locs=nlocs;
                    msoldiff[round]=trimmedMean(locs,0)-optloc; mobjdiff[round]=meanScr(tfunc,locs)-optscr;
                }
                if (data.size()<=10) {
                    out.println("msoldiff="+Arrays.toString(msoldiff));
                    out.println("mobjdiff="+Arrays.toString(mobjdiff));
                    out.println("locs_fin="+Arrays.toString(locs));
                }
                out.println("msoldiff_fin="+msoldiff[ITER]);
                out.println("mobjdiff_fin="+mobjdiff[ITER]);
                out.println("locs_fin="+Arrays.toString(locs));
            }
        }
        out.close();
        System.out.println("time="+(System.currentTimeMillis()-st));
    }
}
