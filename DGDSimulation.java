import java.util.*;
import java.io.*;
public class DGDSimulation {
    private static final String IOFOLDER="..\\..\\DGD_IO\\";
    private static String sign;
    private static double val(double[] p, double x) {
        return p[0]+p[1]*x+p[2]*x*x;
    }
    private static double grad(double[] p, double x) {
        return p[1]+p[2]*2*x;
    }
    private static int[][] graph(String fin) throws IOException {
        BufferedReader in=new BufferedReader(new FileReader(fin));
        sign=in.readLine();
        StringTokenizer tok=new StringTokenizer(sign.split("=")[1].split(",")[0]);
        int N=Integer.parseInt(tok.nextToken());
        int[][] out=new int[N][];
        for (int i=0; i<N; i++) {
            String[] l=in.readLine().split(":");
            int v=Integer.parseInt(l[0]);
            String[] adjs=l[1].split(",");
            out[v]=new int[adjs.length];
            for (int j=0; j<adjs.length; j++)
                out[v][j]=Integer.parseInt(adjs[j]);
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
        int A=0, ITER=10000;
        int N;
        int[][] G; {
        int[][] tG=graph(IOFOLDER+"input_graph.txt");
            N=tG.length;
            G=new int[N+A][];
            System.arraycopy(tG,0,G,0,N);
            int[] tmp=new int[N];
            for (int i=0; i<N; i++) tmp[i]=i;
            for (int i=N; i<N+A; i++)
                G[i]=tmp;
        }
        for (int[] r:G) System.out.println(Arrays.toString(r));
        double[][] Ps=new double[N][];
        double[] locs0=new double[N];
        String tsign;
        {
            BufferedReader in=new BufferedReader(new FileReader(IOFOLDER+"input_testcase.txt"));
            tsign=in.readLine();
            for (int i=0; i<N; i++) {
                String[] l=in.readLine().split(":");
                int v=Integer.parseInt(l[0]);
                double[] p=new double[3];
                String[] info=l[1].split(",");
                for (int j=0; j<3; j++)
                    p[j]=Double.parseDouble(info[j]);
                Ps[v]=p;
            }
            String[] info=in.readLine().split(",");
            for (int j=0; j<N; j++)
                locs0[j]=Double.parseDouble(info[j]);
        }
        for (double[] p:Ps)
            System.out.println(Arrays.toString(p));
        System.out.println(Arrays.toString(locs0));
        double[] tP=new double[3];
        for (int i=0; i<N; i++)
            for (int j=0; j<3; j++)
                tP[j]+=Ps[i][j];
        for (int j=0; j<3; j++) tP[j]/=N;
        double optloc=-tP[1]/(2*tP[2]), optscr=val(tP,optloc);
        System.out.println(Arrays.toString(tP)+"\n"+optloc+","+optscr);

        PrintWriter out=new PrintWriter(new FileWriter(IOFOLDER+"stats_out.txt"));
        out.println(sign);
        out.println(tsign);
        out.printf("A,ITER=%d,%d%n",A,ITER);
        out.printf("optloc,optscr=%f,%f%n",optloc,optscr);

        for (double T:new double[] {0.01,0.005,0.002,0.001}) {
            double[] locs=new double[N];
            System.arraycopy(locs0,0,locs,0,N);
            double[] meanloc=new double[ITER+1], meandiff=new double[ITER+1];
            meanloc[0]=trimmedMean(locs,0); meandiff[0]=meanScr(tP,locs);
            for (int round=1; round<=ITER; round++) {
                List<List<Double>> receivedVals=new ArrayList<>();
                for (int i=0; i<N; i++)
                    receivedVals.add(new ArrayList<>());
                for (int i=0; i<N+A; i++)
                    for (int j:G[i])
                        receivedVals.get(j).add(i<N?locs[i]:1000_000.0);
                double[] nlocs=new double[N];
                for (int i=0; i<N; i++) {
                    double[] tmp=new double[receivedVals.get(i).size()];
                    for (int j=0; j<tmp.length; j++)
                        tmp[j]=receivedVals.get(i).get(j);
                    nlocs[i]=trimmedMean(tmp,A)-T*grad(Ps[i],locs[i]);
                }
                locs=nlocs;
                meanloc[round]=trimmedMean(locs,0); meandiff[round]=meanScr(tP,locs)-optscr;
            }
            out.println("T="+T);
            StringBuilder tmp=new StringBuilder("meanloc=");
            for (int i=0; i<=ITER; i++)
                tmp.append(meanloc[i]).append(i<ITER?",":"\n");
            out.print(tmp);
            tmp=new StringBuilder("meandiff=");
            for (int i=0; i<=ITER; i++)
                tmp.append(meandiff[i]).append(i<ITER?",":"\n");
            out.print(tmp);
        }
        out.close();
        System.out.println("time="+(System.currentTimeMillis()-st));
    }
}
