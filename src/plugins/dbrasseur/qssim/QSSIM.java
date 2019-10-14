package plugins.dbrasseur.qssim;

import icy.image.IcyBufferedImage;
import plugins.dbrasseur.qssim.quaternion.Quat;

public class QSSIM {
    public static double computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg){return computeQSSIM(ref, deg, 0.01, 0.03);}
    public static double computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2){return computeQSSIM(ref, deg, K1, K2, 1.5);}
    public static double computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double std){return computeQSSIM(ref, deg, K1, K2, std, std, 11);}
    public static double computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double stdX, double stdY, int kSize)
    {
        int w = ref.getWidth();
        int h = ref.getHeight();
        //Images checkup
        //TODO

        //Convert image to Quaternions
        Quat[] refq = imageToQuaternionf(ref);
        Quat[] degq = imageToQuaternionf(deg);

        Quat C1=new Quat(K1*K1, 0, 0, 0);
        Quat C2=new Quat(K2*K2, 0, 0, 0);

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
