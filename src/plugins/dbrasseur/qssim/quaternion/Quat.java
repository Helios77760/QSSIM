package plugins.dbrasseur.qssim.quaternion;

public class Quat {
    private double re;
    private Vec3 im;

    public Quat(double a, double b, double c, double d)
    {
        this(a, new Vec3(b,c,d));
    }

    public Quat(double re, Vec3 im) {
        this.re = re;
        this.im = im;
    }

    public Quat(){
        this(0.0, new Vec3());
    }

    public Quat(Quat a)
    {
        this(a.re, new Vec3(a.im));
    }

    public double getRe() {
        return re;
    }

    public void setRe(double re) {
        this.re = re;
    }

    public Vec3 getIm() {
        return im;

    }

    public Quat add(Quat o)
    {
        re += o.re;
        im.add(o.im);
        return this;
    }

    public static Quat add(Quat a, Quat b)
    {
        Quat res = new Quat(a);
        res.add(b);
        return res;
    }

    public Quat mul(double a)
    {
        re*=a;
        im.mul(a);
        return this;
    }

    public static Quat mul(Quat a, double b)
    {
        Quat res = new Quat(a);
        res.mul(b);
        return res;
    }

    public static Quat mul(double a, Quat b)
    {
        return mul(b, a);
    }

    public static Quat mul(Quat a, Quat b)
    {
        return new Quat(a.re*b.re - a.im.dot(b.im), Vec3.mul(a.re, b.im).add(Vec3.mul(b.re, a.im)).add(Vec3.cross(a.im, b.im)));
    }

    public Quat sub(Quat o)
    {
        re-=o.re;
        im.sub(o.im);
        return this;
    }

    public static Quat sub(Quat a, Quat b)
    {
        Quat res = new Quat(a);
        res.sub(b);
        return res;
    }

    public static double dot(Quat a, Quat b)
    {
        return a.dot(b);
    }

    public double dot(Quat o)
    {
        return re*o.re+im.dot(o.im);
    }

    public Quat inverse()
    {
        return (conjugate().mul(1.0/squaredNorm()));
    }

    public Quat conjugate()
    {
        return new Quat(re, -im.x, -im.y, -im.z);
    }

    public double norm()
    {
        return Math.sqrt(squaredNorm());
    }

    public double squaredNorm(){return re*re+im.x*im.x+im.y*im.y+im.z*im.z;}


}