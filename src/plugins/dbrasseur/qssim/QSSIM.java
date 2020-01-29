package plugins.dbrasseur.qssim;

import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.type.DataType;
import plugins.dbrasseur.qssim.quaternion.Quat;

//Main implementation class
public class QSSIM {
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, int[] map_size, boolean downsample)throws IllegalArgumentException{return computeQSSIM(ref, deg, 0.01, 0.03, map_size, downsample);}
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, int[] map_size, boolean downsample)throws IllegalArgumentException{return computeQSSIM(ref, deg, K1, K2, 255, map_size, downsample);}
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double dynamicRange, int[] map_size, boolean downsample)throws IllegalArgumentException{return computeQSSIM(ref, deg, K1, K2, dynamicRange, 1.5, map_size, downsample);}
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double dynamicRange, double std, int[] map_size, boolean downsample)throws IllegalArgumentException{return computeQSSIM(ref, deg, K1, K2, dynamicRange, std, std, map_size, downsample);}
    
    /**
     * Computes the QSSIM map between 2 images
     * @param ref Reference image
     * @param deg Degraded image
     * @param K1 Regularisation constant
     * @param K2 Regularisation constant
     * @param dynamicRange dynamic range of the input image
     * @param stdX Standard deviation in the X direction of the gaussian kernel window
     * @param stdY Standard deviation in the Y direction of the gaussian kernel window
     * @param map_size Output for the final map size (potentially downsampled)
     * @return QSSIM map array
     * @throws IllegalArgumentException if the images are not proper
     */
    public static double[] computeQSSIM(IcyBufferedImage ref, IcyBufferedImage deg, double K1, double K2, double dynamicRange, double stdX, double stdY, int[] map_size, boolean downsample) throws IllegalArgumentException
    {
        //Images checkup
        if(ref == null || deg == null)
        {
            throw new IllegalArgumentException("Null image");
        }else{
            if(ref.getSizeC() != deg.getSizeC() || ref.getSizeX() != deg.getSizeX() || ref.getSizeY() != deg.getSizeY())
            {
                throw new IllegalArgumentException("Images dimensions don't match");
            }else if(ref.getSizeX() == 0 || ref.getSizeY() == 0){
                throw new IllegalArgumentException("Bad images size : X:"+ref.getSizeX()+" Y:"+ref.getSizeY()+" C:"+ref.getSizeC());
            }else if(ref.getDataType_() != deg.getDataType_()){
                throw new IllegalArgumentException("Mismatched image Data Types");
            }else if(K1 <= 0 || K2 <= 0)
            {
                throw new IllegalArgumentException("Invalid regularisation parameter, must be >= 0");
            }else if(dynamicRange < 1)
            {
                throw new IllegalArgumentException("Invalid dynamic range, must be >= 1, rapported to an 8 bits integer");
            }else if(stdX < 0 || stdY < 0)
            {
                throw new IllegalArgumentException("Invalid standard deviation for the kernel, must be >= 0");
            }else if(map_size == null)
            {
            	throw new IllegalArgumentException("Map size output is null");
            }
        }
        
        double L = dynamicRange;
        //Get the datatype to adjust the range
        DataType imgType = ref.getDataType_();
        switch(imgType){
            case UBYTE:
            case BYTE:
                L/=0xFF;
                break;
            case USHORT:
            case SHORT:
                L/=0xFFFF;
                break;
            case UINT:
            case INT:
                L/=0xFFFFFFFF;
                break;
            case ULONG:
            case LONG:
                L/=0xFFFFFFFFFFFFFFFFL;
                break;
            default:
                //Nothing to do
                break;
        }
        
        //Get number of processors for multi-threading
        final int nbProc = Runtime.getRuntime().availableProcessors() <= 0 ? 1 : Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[nbProc];

        //Define regularisation constants
        final Quat C1=new Quat(K1*K1*L*L, 0, 0, 0);
        final Quat C2=new Quat(K2*K2*L*L, 0, 0, 0);

        //Convert image to Quaternions
        Quat[] refq = imageToQuaternionf(IcyBufferedImageUtil.convertToType(ref,DataType.DOUBLE, true), map_size, downsample);
        Quat[] degq = imageToQuaternionf(IcyBufferedImageUtil.convertToType(deg,DataType.DOUBLE, true), map_size, downsample);
        
        //Get the weight kernel
        final double[][] kernel = createGaussian(stdX, stdY);
        final int padding = kernel.length/2;
        
        //Modify the sizes to accomodate the window
        int imw = map_size[0];
        int w = map_size[0] =  (map_size[0]-2*padding);
        int h = map_size[1] =  (map_size[1]-2*padding);
        
        //Compute all the windows in parallel
        double[] qssim_map = new double[w*h];
        for(int p=0; p<nbProc; p++)
        {
            final int currP=p;
            threads[p] = new Thread(()->{
                for(int i=currP*w*h/nbProc; i<((currP+1)*w*h)/nbProc; i++)
                {
                    int x = i/w;
                    int y = i%w;

                    Quat[][] patchref = extract2D(refq, x+padding, y+padding, padding, imw);
                    Quat[][] patchdeg = extract2D(degq, x+padding, y+padding, padding, imw);
                    qssim_map[x*w+y] = computeQSSIMPatch(patchref, patchdeg, C1, C2, kernel);
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
    }

    /**
     * Extracts the window at the given center into a 2D array
     * @param im Flat image
     * @param centerX X coordinate of the center of the window 
     * @param centerY Y coordinate of the center of the window 
     * @param radius Radius of the window
     * @param w Width of the flatten image
     * @return 2D window
     */
    private static Quat[][] extract2D(Quat[] im, int centerX, int centerY, int radius, int w)
    {
        Quat[][] res = new Quat[radius*2+1][radius*2+1];
        for(int x=-radius; x <= radius; x++)
        {
        	System.arraycopy(im, (centerX+x)*w+centerY-radius, res[x+radius], 0, radius*2+1);
        }
        return res;
    }

    /**
     * Creates a gaussian window
     * @param stdX Standard deviation in the X direction
     * @param stdY Standard deviation in the Y direction
     * @return 2D Gaussian window
     */
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
            	kernel[i][j] = rawY[i]*rawX[j];
                sum+= kernel[i][j];
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

    /**
     * Center a 1D array
     * @param arr Input array
     * @param size new size for the array
     * @return 0 padded array of the input array
     */
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

    /**
     * Creates a Gaussian 1D array
     * @param std Standard deviation
     * @return 1D Gaussian array
     */
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
    
    /**
     * Computes the QSSIM value of a given patch
     * @param ref Reference patch
     * @param deg Degraded patch
     * @param C1 Regularisation constant
     * @param C2 Regularisation constant
     * @param weight Gaussian window
     * @return QSSIM value
     */
    private static double computeQSSIMPatch(Quat[][] ref, Quat[][] deg, Quat C1, Quat C2, double[][] weight)
    {
        //Compute mean luminance : muq (13)
        Quat muq_ref = meanLuminance(ref, weight);
        Quat muq_deg = meanLuminance(deg, weight);
        
        //Substract the mean luminance from the image : acq (15)
        Quat[][] acq_ref = chrominance(ref, muq_ref);
        Quat[][] acq_deg = chrominance(deg, muq_deg);

        //Compute color contrast sigmaq(17)
        Quat sigmaq_ref_sq = correlation(acq_ref, acq_ref, weight);
        Quat sigmaq_deg_sq = correlation(acq_deg, acq_deg, weight);

        //Compute cross correlation (20)
        Quat sigmaq_ref_deg = correlation(acq_ref, acq_deg, weight);

        //Compute QSSIM (between 20 and 21)
        Quat num1 = Quat.mul(2.0, Quat.mul(muq_ref, muq_deg.conjugate())).add(C1);
        Quat num2 = Quat.mul(2.0,sigmaq_ref_deg).add(C2);
        Quat den1 = Quat.mul(muq_ref, muq_ref.conjugate()).add(Quat.mul(muq_deg, muq_deg.conjugate())).add(C1);
        Quat den2 = Quat.add(sigmaq_deg_sq, sigmaq_ref_sq).add(C2);

        return Quat.mul(Quat.mul(num1,num2), Quat.mul(den1, den2).inverse()).norm();
    }

    /**
     * Correlation between chrominances
     * @param win1 First window
     * @param win2 Second window
     * @param weight Gaussian Window
     * @return Correlation
     */
    private static Quat correlation(Quat[][] win1, Quat[][] win2, double[][] weight) {
        Quat result = new Quat();
        for(int i=0; i<weight.length; i++)
        {
            for(int j=0; j<weight[i].length; j++)
            {
                result.add(Quat.mul(win1[i][j], win2[i][j].conjugate()).mul(weight[i][j]));
            }
        }
        return result;
    }

    /**
     * Computes the chrominance
     * @param im window
     * @param mean luminance
     * @return chrominance window
     */
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

    /**
     * Computes the mean luminance
     * @param im window
     * @param weight gaussian window
     * @return mean luminance
     */
    private static Quat meanLuminance(Quat[][] im, double[][] weight) {
        Quat sum= new Quat();
        for(int i=0; i<im.length; i++)
        {
            for(int j=0; j<im[i].length; j++)
            {
                sum.add(Quat.mul(im[i][j], weight[i][j]));
            }
        }
        return sum;
    }

    /**
     * Converts an image into its quaternion array form
     * @param img Image
     * @param newSizes Output the new size in case of downsampling
     * @param autoDownsample Enables automatic downsampling
     * @return Quaternion version of the image
     */
    public static Quat[] imageToQuaternionf(IcyBufferedImage img, int[] newSizes, boolean autoDownsample)
    {
        double[][] ref_pixels = img.getDataXYCAsDouble();
        
        int w=img.getWidth();
        int h=img.getHeight();
        double[][] base = new double[3][w*h];
        
        //Automatic Downsampling
        int f = 1;
        if(autoDownsample)
        {
        	f = (int) Math.max(1, Math.round(Math.min(w, h)/256.0));
            double filtervalue = 1.0/(f*f);
            int foff = f%2 == 0 ? -f/2+1 : -f/2;
            //Mirror filtering
            for(int c=0; c<3; c++)
            {
                for(int i=0; i<img.getHeight(); i++)
                {
                    for(int j=0; j<img.getWidth(); j++)
                    {
                    	base[c][i*w+j] = 0;
                        for(int fi=0; fi<f; fi++)
                        {
                            int fx = i+foff+fi;
                            if(fx < 0)
                            {
                                fx = -fx-1;
                            }else if(fx >= h)
                            {
                            	fx = 2*h-fx-1;
                            }
                            for(int fj=0; fj<f; fj++)
                            {
                            	int fy = j+foff+fj;
                                if(fy < 0)
                                {
                                    fy = -fy-1;
                                }else if(fy >= w)
                                {
                                	fy = 2*w-fy-1;
                                }
                                base[c][i*w+j] += ref_pixels[c][fx*w+fy]*filtervalue;
                            }
                        }
                    }
                }
            }
        }else
        {
        	base = ref_pixels;
        }
       
        
        newSizes[0] = w/f;
        newSizes[1] = h/f;
        Quat[] refq = new Quat[newSizes[1]*newSizes[0]];
        int out=0;
        for(int i=0; i<newSizes[1]; i++)
        {
        	for(int j=0; j<newSizes[0]; j++)
        	{
                int index = i * f * w + j * f;
                refq[out++] = new Quat(0.0, base[0][index]/Math.sqrt(3), base[1][index]/Math.sqrt(3), base[2][index]/Math.sqrt(3));
        	}
        }
        
        return refq;
    }


    /**
     * Computes the mean of the array
     * @param arr Array
     * @return Mean
     */
    public static double meanArray(double[] arr)
    {
        if(arr.length == 0) return 0;
        double sum=0;
        for(double d : arr) sum+=d;
        return sum/arr.length;
    }
}
