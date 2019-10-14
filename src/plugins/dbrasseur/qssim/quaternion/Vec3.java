package plugins.dbrasseur.qssim.quaternion;

public class Vec3 {
    public double x;
    public double y;
    public double z;

    public Vec3(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vec3()
    {
        this(0.0, 0.0, 0.0);
    }

    public Vec3(Vec3 o)
    {
        this(o.x, o.y, o.z);
    }

    public Vec3 add(Vec3 o)
    {
        x+=o.x;
        y+=o.y;
        z+=o.z;
        return this;
    }

    public Vec3 sub(Vec3 o)
    {
        x-=o.x;
        y-=o.y;
        z-=o.z;
        return this;
    }

    public Vec3 mul(double a)
    {
        x*=a;
        y*=a;
        z*=a;
        return this;
    }

    public static Vec3 add(Vec3 a, Vec3 b)
    {
        Vec3 res = new Vec3(a);
        return res.add(b);
    }

    public static Vec3 sub(Vec3 a, Vec3 b)
    {
        Vec3 res = new Vec3(a);
        return res.sub(b);
    }

    public static Vec3 mul(double a, Vec3 v)
    {
        Vec3 res = new Vec3(v);
        return res.mul(a);
    }

    public double dot(Vec3 o)
    {
        return x*o.x+y*o.y+z*o.z;
    }

    public static double dot(Vec3 a, Vec3 b)
    {
        return a.dot(b);
    }

    public static Vec3 cross(Vec3 a, Vec3 b)
    {
        return new Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
    }
}
