package plugins.dbrasseur.qssim;

import icy.gui.dialog.MessageDialog;
import icy.gui.frame.GenericFrame;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.sequence.Sequence;
import icy.type.DataType;
import plugins.adufour.ezplug.*;

import javax.swing.*;
import java.awt.*;

/**
 * Implementation of "Quaternion Structural Similarity: A New Quality Index for Color Images" by Kolaman and Yadid
 * @author Dylan Brasseur
 * @version 1.3
 *
 */
public class QSSIMPlug extends EzPlug{

	private EzVarSequence	EzSrcSeq;				//Source Sequence
    private EzVarSequence   EzDegSeq;               //Degraded Sequence
    private EzVarInteger	EzDivider;				//Number of degraded images per reference
	private EzVarDouble		EzL;					//Dynamic range parameter
	private EzVarDouble		EzK1;					//K1 parameter
	private EzVarDouble		EzK2;					//K2 parameter
	private EzVarDouble		EzSigmaX;				//SigmaX parameter
	private EzVarDouble		EzSigmaY;				//SigmaY parameter
	private EzVarBoolean	EzDownsample;			//Enables downsampling


	@Override
	public void clean() {

	}

	@Override
	protected void execute() {
        Sequence src = EzSrcSeq.getValue();
        Sequence deg = EzDegSeq.getValue();
        Integer divider = EzDivider.getValue();
        if(divider == 0)
        {
        	//auto split
        	divider = Math.max(1, deg.getNumImage()/src.getNumImage());
        }
        if(src == null || deg==null)
		{
			MessageDialog.showDialog("Please open an image first.", MessageDialog.ERROR_MESSAGE);
			return;
		}
        double K1 = EzK1.getValue();
        double K2 = EzK2.getValue();
        double stdX = EzSigmaX.getValue();
        double stdY = EzSigmaY.getValue();
        double L = EzL.getValue();
        boolean downsample = EzDownsample.getValue();
        Sequence res = new Sequence();
        for(int srci=0; srci < src.getNumImage(); srci++)
        {
        	IcyBufferedImage srcImg = src.getImage(srci, 0);
        	boolean grey = srcImg.getSizeC() < 3;
    		if(grey) srcImg = toGreyscale(srcImg);
    		
    		for(int i=0; i<divider; i++)
    		{
    			IcyBufferedImage degImg = deg.getImage(i, srci);
    			if(srcImg.getWidth() != degImg.getWidth() || srcImg.getHeight() != degImg.getHeight())
    			{
    				MessageDialog.showDialog("Image dimensions ", MessageDialog.ERROR_MESSAGE);
    				return;
    			}
    			
    			if(grey) degImg = toGreyscale(degImg);
    			int[] map_size = new int[2];
    			double[] qssim_map = QSSIM.computeQSSIM(srcImg, degImg, K1, K2, L, stdX, stdY, map_size, downsample);
    			double mqssim = QSSIM.meanArray(qssim_map);
    			new AnnounceFrame("The mean QSSIM between "+ src.getName()  + "_" + srci + " and "+deg.getName() + "_" + i + " is "+mqssim);
    			System.out.println("The mean QSSIM between "+ src.getName()  + "_" + srci + " and "+deg.getName() + "_" + i + " is "+mqssim);
    			IcyBufferedImage map = new IcyBufferedImage(map_size[0], map_size[1], 1, DataType.DOUBLE);
    			map.setDataXY(0, qssim_map);
    			map = IcyBufferedImageUtil.convertToType(map, DataType.UBYTE, true);
    			res.addImage(map);
    			
    		}
        }
		
		res.setName("QSSIM : "+src.getName()+" / " + deg.getName());
		addSequence(res);
	}

	private IcyBufferedImage toGreyscale(IcyBufferedImage srcImg) {
		System.out.println("Converting image from greyscale as it doesn't have at least 3 components");
		IcyBufferedImage newSrc = new IcyBufferedImage(srcImg.getWidth(), srcImg.getHeight(), 3, DataType.DOUBLE);
		newSrc.setDataXY(0, srcImg.getDataXY(0));
		newSrc.setDataXY(1, srcImg.getDataXY(0));
		newSrc.setDataXY(2, srcImg.getDataXY(0));
		return newSrc;
	}

	@Override
	protected void initialize() {
		super.setTimeDisplay(true);
		
		EzSrcSeq = new EzVarSequence("Reference");
		EzSrcSeq.setToolTipText("Reference sequence");
		
		EzDegSeq = new EzVarSequence("Degraded");
		EzDegSeq.setToolTipText("Degraded sequence");
		
		EzDivider = new EzVarInteger("nDeg/ref", 0, Integer.MAX_VALUE, 1);
		EzDivider.setToolTipText("Number of degraded images per reference image, 0 = auto");
		
		EzL = new EzVarDouble("Dynamic Range",255, 1, Integer.MAX_VALUE, 1);
		EzL.setToolTipText("Dynamic range of the image (typically 255 for 8 bits/pixel)");
		
		EzK1 = new EzVarDouble("K1", 0.01, Double.MIN_NORMAL, Double.MAX_VALUE, 0.0001);
		EzK1.setToolTipText("Regularization parameter, should be as small as possible but > 0, defaults to 0.01");
		
		EzK2 = new EzVarDouble("K2", 0.03, Double.MIN_NORMAL, Double.MAX_VALUE, 0.0001);
		EzK2.setToolTipText("Regularization parameter, should be as small as possible but > 0, defaults to 0.03");
		
		EzSigmaX = new EzVarDouble("\u03C3x", 1.5, 0, 1024, 0.1);
		EzSigmaX.setToolTipText("Standard deviation on X of the kernel used in the window, defaults to 1.5");
		
		EzSigmaY = new EzVarDouble("\u03C3y", 1.5, 0, 1024, 0.1);
		EzSigmaY.setToolTipText("Standard deviation on Y of the kernel used in the window, defaults to 1.5");
		
		EzDownsample = new EzVarBoolean("Downsample", true);
		EzDownsample.setToolTipText("Automatically downsample the image, defaults to true");
		
		EzGroup parameters = new EzGroup("Parameters", EzL, EzK1, EzK2, EzSigmaX, EzSigmaY, EzDownsample);
		super.addEzComponent(EzSrcSeq);
		super.addEzComponent(EzDegSeq);
		super.addEzComponent(EzDivider);
		super.addEzComponent(parameters);
		EzButton ezHelp = new EzButton("Informations about QSSIM", (e -> printInformations()));
		super.addEzComponent(ezHelp);
	}

	@SuppressWarnings("deprecation")
	private void printInformations() {
		String    title   = "QSSIM informations";
		JTextPane message = new JTextPane();
		message.setEditable(false);
		message.setContentType("text/html");
		message.setText(
				"<p>" +
						"The QSSIM (Quaternion Structural Similarity) is an index measuring the structural similarity between two color images. " +
						"It is valued between 0 and 1. When two images are nearly identical, " +
						"their QSSIM is close to 1." +
						"</p>" +
						"<p>" +
						"Formula for computing the QSSIM between two sequences seq1 and seq2 at a given pixel " +
						"or voxel <tt>P</tt>:" +
						"<pre>" +
						"            |   2*mu1(P)*mu2(P) + C1         2*\u03C312(P) + C2     |\n" +
						"  QSSIM(P) =| ------------------------ x ----------------------|\n" +
						"            | mu1(P)^2 + mu2(P)^2 + C1   \u03C31(P)^2 + \u03C32(P)^2 + C2|\n" +
						"</pre>" +
						"With:" +
						"<ul>" +
						"<li><tt>mu1(P)</tt> and <tt>mu2(P)</tt>: mean value of seq1 and seq2 computed " +
						"over a small XY window located around <tt>P</tt></li>" +
						"<li><tt>\u03C31(P)</tt> and <tt>\u03C32(P)</tt>: standard deviation of seq1 and seq2 computed " +
						"over the same window</li>" +
						"<li><tt>\u03C312(P)</tt>: color cross correlation between seq1 and seq2 computed over the same window</li>" +
						"<li><tt>C1 = (K1*L)^2</tt>: regularization constant (should be as small as possible)</li>" +
						"<li><tt>C2 = (K2*L)^2</tt>: regularization constant (should be as small as possible)</li>" +
						"<li><tt>K1</tt>, <tt>K2</tt>: regularization parameters (must be >0)</li>" +
						"<li><tt>L</tt>: dynamic range of the pixel values (example: <tt>L=255</tt> " +
						"if the sequence is 8 bit encoded)</li>" +
						"</ul>" +
						"</p>" +
						"<p>" +
						"The default window is a Gaussian window with standard deviation 1.5 along both " +
						"the X and the Y axis." +
						"</p>" +
						"<p>" +
						"Reference:" +
						"<quote>" +
						" A. Kolaman and O. Yadid-Pecht, \u00AB Quaternion Structural Similarity: A New Quality Index for Color Images \u00BB, IEEE Trans. on Image Process., vol. 21, n\u00B0 4, p. 1526-1536, avr. 2012." +
						"</quote>" +
						"</p>" +
						"<p>" +
						"This current implementation sticks as much as possible to the Matlab QSSIM " +
						"implementation provided by these authors at:<br/>" +
						"<a href=\"http://www.ee.bgu.ac.il/~kolaman/QSSIM\">" +
						"http://www.ee.bgu.ac.il/~kolaman/QSSIM" +
						"</a><br/>This implementation is also inspired from the Icy SSIM Toolbox Plugin by Yoann Le Montagner." +
						"</p>"
		);
		Dimension dim = message.getPreferredSize();
		dim.setSize(600, dim.getHeight()+100);
		message.setPreferredSize(dim);
		JScrollPane scroll = new JScrollPane(message);
		dim = scroll.getPreferredSize();
		dim.setSize(600, 500);
		scroll.setPreferredSize(dim);
		GenericFrame infoFrame = new GenericFrame(title, scroll);
		infoFrame.addToMainDesktopPane();
		infoFrame.setVisible(true);
		infoFrame.requestFocus();
	}
}