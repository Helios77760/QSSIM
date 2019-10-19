package plugins.dbrasseur.qssim;

import icy.image.IcyBufferedImage;
import plugins.dbrasseur.qssim.quaternion.Quat;

public class QSSIM {
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg){return computeQSSIM(ref, deg, 0.01, 0.03);}
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2){return computeQSSIM(ref, deg, K1, K2, 1.5);}
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double std){return computeQSSIM(ref, deg, K1, K2, std, std);}
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double stdX, double stdY)
    {
        int w = ref.getWidth();
        int h = ref.getHeight();
        //Images checkup
        //TODO

        //Convert image to Quaternions
        Quat[] refq = imageToQuaternionf(ref);
        Quat[] degq = imageToQuaternionf(deg);

        final Quat C1=new Quat(K1*K1, 0, 0, 0);
        final Quat C2=new Quat(K2*K2, 0, 0, 0);

        double[] qssim_map = new double[w*h];
        final double[][] kernel = createGaussian(stdX, stdY);
        final Quat[][] paddedRef = pad_mirror(refq, w, h, kernel.length);
        final Quat[][] paddedDeg = pad_mirror(degq, w, h, kernel.length);
        final int nbProc = Runtime.getRuntime().availableProcessors() <= 0 ? 1 : Runtime.getRuntime().availableProcessors();
        int padding = kernel.length/2;
        Thread[] threads = new Thread[nbProc];
        for(int p=0; p<nbProc; p++)
        {
            final int currP=p;
            threads[p] = new Thread(()->{
                for(int i=(currP*w*h)/nbProc; i<((currP+1)*w*h)/nbProc; i++)
                {
                    int x = i/w;
                    int y = i%w;

                    Quat[][] patchref = extract2D(paddedRef, x+padding, y+padding, padding);
                    Quat[][] patchdeg = extract2D(paddedDeg, x+padding, y+padding, padding);
                    qssim_map[x*w+y] = computeQSSIMPatch(patchref, patchdeg, kernel.length, kernel.length, C1, C2, kernel);
                }
            });
            threads[p].start();
        }
        for(Thread T : threads)
        {
            try {
                T.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        return qssim_map;
/*

        //Compute mean luminance : muq (13)
        Quat muq_ref = meanLuminancef(refq);
        Quat muq_deg = meanLuminancef(degq);

        //Substract the mean luminance from the image : acq (15)
        Quat[] acq_ref = chrominancef(refq, muq_ref);
        Quat[] acq_deg = chrominancef(degq, muq_deg);

        //Compute color contrast sigmaq(17)
        Quat sigmaq_ref_sq = contrastf(acq_ref, w, h);
        Quat sigmaq_deg_sq = contrastf(acq_deg, w, h);

        //Compute cross correlation (20)
        Quat sigmaq_ref_deg = crossCorrelationf(acq_ref, acq_deg, w, h);

        //Compute QSSIM (between 20 and 21)
        Quat num1 = Quat.mul(2.0, Quat.mul(muq_ref, muq_deg.conjugate())).add(C1);
        Quat num2 = Quat.mul(2.0,sigmaq_ref_deg).add(C2);
        Quat den1 = Quat.mul(muq_ref, muq_ref.conjugate()).add(Quat.mul(muq_deg, muq_deg.conjugate())).add(C1);
        Quat den2 = Quat.add(sigmaq_deg_sq, sigmaq_ref_sq).add(C2);

        return Quat.mul(Quat.mul(num1,num2), Quat.mul(den1, den2).inverse()).norm();
*/
    }

    private static Quat[][] extract2D(Quat[][] im, int centerX, int centerY, int radius)
    {
        Quat[][] res = new Quat[radius*2+1][radius*2+1];
        for(int x=-radius; x <= radius; x++)
        {
            for(int y=-radius; y <= radius; y++)
            {
                res[x+radius][y+radius] = im[centerX+x][centerY+y];
            }
        }
        return res;
    }

    private static Quat[][] pad_mirror(Quat[] im, int w, int h, int length) {
           if(length<=1)
           {
               Quat[][] res = new Quat[h][w];
               for(int i=0; i<h; i++)
               {
                   if (w >= 0) System.arraycopy(im, i * w, res[i], 0, w);
               }
               return res;
           }else
           {
               assert(length%2 == 1);
               int padding = length/2;
               Quat[][] res = new Quat[h+2*padding][w+2*padding];
               for(int i=0; i<h; i++)
               {
                   if (w >= 0) System.arraycopy(im, i * w, res[i+padding], padding, w);
                   for(int j=0; j<padding; j++)
                   {
                       res[i+padding][j] = new Quat(im[i*w+(padding-j-1)]);
                       res[i+padding][padding+w+j] = new Quat(im[i*w+w-j-1]);
                   }
               }
               for(int i=0; i< padding; i++)
               {
                   System.arraycopy(res[2*padding-1-i], 0, res[i], 0, w+2*padding);
                   System.arraycopy(res[h+padding-i-1], 0, res[i+padding+h], 0, w+2*padding);
               }
               return res;
           }
    }

    private static double[][] createGaussian(double stdX, double stdY) {
        double[] rawX = createGaussian1D(stdX);
        double[] rawY = createGaussian1D(stdY);
        int size = rawX.length > rawY.length ? rawX.length : rawY.length;
        double[][] kernel = new double[size][size];


        if(size > rawX.length){ // Center 1D filter
            rawX = center1D(rawX, size);
        }else if(size > rawY.length){
            rawY = center1D(rawY, size);
        }
        double sum=0;
        for(int i=0; i<rawY.length; i++){
            for(int j=0; j<rawX.length; j++)
            {
                sum+= kernel[i][j] = rawY[i]*rawX[i];
            }
        }
        for(int i=0; i<rawY.length; i++){
            for(int j=0; j<rawX.length; j++)
            {
                kernel[i][j]/=sum;
            }
        }
        return kernel;
    }

    private static double[] center1D(double[] arr, int size) {
        int pad = (size-arr.length)/2;
        double[] newRawX = new double[size];
        for(int i=0; i<size; i++){
            if(i-pad < 0 || i-pad >= arr.length)
            {
                newRawX[i]=0;
            }else
            {
                newRawX[i]=arr[i-pad];
            }
        }
        return newRawX;
    }

    private static double[] createGaussian1D(double std)
    {
        if (std < 1.0e-10)
        {
            return new double[] { 1 };
        }

        double sigma2 = std * std;
        int k = (int) Math.ceil(std * 3.0f);

        int width = 2 * k + 1;

        double[] kernel = new double[width];
        double sum=0;
        for (int i = -k; i <= k; i++){
            sum += kernel[i + k] = 1.0 / (Math.sqrt(2 * Math.PI) * std * Math.exp(((i * i) / sigma2) * 0.5f));
        }
        for(int i=0; i<kernel.length; i++)
            kernel[i]/=sum;
        return kernel;
    }

    private static double computeQSSIMPatch(Quat[][] ref, Quat[][] deg, int w, int h, Quat C1, Quat C2, double[][] weight)
    {
        //Compute mean luminance : muq (13)
        Quat muq_ref = meanLuminance(ref, weight);
        Quat muq_deg = meanLuminance(deg, weight);

        //Substract the mean luminance from the image : acq (15)
        Quat[][] acq_ref = chrominance(ref, muq_ref);
        Quat[][] acq_deg = chrominance(deg, muq_deg);

        //Compute color contrast sigmaq(17)
        Quat sigmaq_ref_sq = contrast(acq_ref, w, h, weight);
        Quat sigmaq_deg_sq = contrast(acq_deg, w, h, weight);

        //Compute cross correlation (20)
        Quat sigmaq_ref_deg = crossCorrelation(acq_ref, acq_deg, w, h, weight);

        //Compute QSSIM (between 20 and 21)
        Quat num1 = Quat.mul(2.0, Quat.mul(muq_ref, muq_deg.conjugate())).add(C1);
        Quat num2 = Quat.mul(2.0,sigmaq_ref_deg).add(C2);
        Quat den1 = Quat.mul(muq_ref, muq_ref.conjugate()).add(Quat.mul(muq_deg, muq_deg.conjugate())).add(C1);
        Quat den2 = Quat.add(sigmaq_deg_sq, sigmaq_ref_sq).add(C2);

        return Quat.mul(Quat.mul(num1,num2), Quat.mul(den1, den2).inverse()).norm();
    }

    private static Quat crossCorrelation(Quat[][] ref, Quat[][] deg, int w, int h, double[][] weight) {
        Quat result = new Quat();
        for(int i=0; i<ref.length; i++)
        {
            for(int j=0; j<ref[i].length; j++)
            {
                result.add(Quat.mul(ref[i][j], deg[i][j].conjugate()).mul(weight[i][j]));
            }

        }
        return result.mul(1.0/((w-1)*(h-1)));
    }

    private static Quat contrast(Quat[][] ac, int w, int h, double[][] weight) {
        Quat sum=new Quat();
        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++)
            {
                sum.add(Quat.mul(ac[i][j], ac[i][j].conjugate()).mul(weight[i][j]));
            }
        }
        return sum.mul(1.0/((w-1)*(h-1)));
    }

    private static Quat[][] chrominance(Quat[][] im, Quat mean) {
        Quat[][] res = new Quat[im.length][im[0].length];
        for(int i=0; i<im.length; i++)
        {
            for(int j=0; j<im[0].length; j++)
            {
                res[i][j] = Quat.sub(im[i][j], mean);
            }
        }
        return res;
    }

    private static Quat meanLuminance(Quat[][] im, double[][] weight) {
        Quat sum= new Quat();
        for(int i=0; i<im.length; i++)
        {
            for(int j=0; j<im[i].length; j++)
            {
                sum.add(Quat.mul(im[i][j], weight[i][j]));
            }
        }
        return sum.mul(1.0/(im.length*im[0].length));
    }

    private static Quat crossCorrelationf(Quat[] acq_ref, Quat[] acq_deg, int w, int h) {
            Quat result = new Quat();
            for(int i=0; i<acq_ref.length; i++)
            {
                result.add(Quat.mul(acq_ref[i], acq_deg[i].conjugate()));
            }
            return result.mul(1.0/((w-1)*(h-1)));
    }

    private static Quat contrastf(Quat[] CV, int w, int h) {
        Quat sigma = new Quat();
        for(Quat q : CV)
        {
            sigma.add(Quat.mul(q, q.conjugate()));
        }
        return sigma.mul(1.0/((w-1)*(h-1)));
    }

    private static Quat[] chrominancef(Quat[] img, Quat mean) {
        Quat[] CV = new Quat[img.length];
        for(int i=0; i < img.length; i++)
        {
            CV[i] = Quat.sub(img[i], mean);
        }
        return CV;
    }

    public static Quat[] imageToQuaternionf(IcyBufferedImage img)
    {
        double[][] ref_pixels = img.getDataXYCAsDouble();
        Quat[] refq = new Quat[img.getWidth()*img.getHeight()];
        for(int i=0; i<img.getWidth()*img.getHeight(); i++)
        {
            refq[i] = new Quat(0.0, ref_pixels[0][i], ref_pixels[1][i], ref_pixels[2][i]);
        }
        return refq;
    }



    private static Quat meanLuminancef(Quat[] img)
    {
        Quat mean = new Quat();
        for(Quat q : img)
        {
            mean.add(q);
        }
        return mean.mul(1.0/img.length);

    }
}
