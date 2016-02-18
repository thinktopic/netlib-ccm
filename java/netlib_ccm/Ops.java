package netlib_ccm;


public class Ops
{
    public static void alphaAXPlusBetaY( int op_len, double alpha,
					 double[] a, int a_offset,
					 double[] x, int x_offset,
					 double beta,
					 double[] y, int y_offset)
    {
	for (int idx = 0; idx < op_len; ++idx) {
	    int y_off = y_offset + idx;
	    y[y_off] = alpha * x[x_offset + idx] * a[a_offset + idx] + beta * y[y_off];
	}
    }

    public static void alphaAY( int op_len, double alpha,
				double[] a, int a_offset,
				double[] y, int y_offset )
    {
	for (int idx = 0; idx < op_len; ++idx) {
	    int y_off = y_offset + idx;
	    y[y_off] = alpha * a[a_offset + idx] * y[y_off];
	}
    }

    public static void axpy( int op_len, double alpha,
			     double[] x, int x_offset,
			     double[] y, int y_offset )
    {
	for (int idx = 0; idx < op_len; ++idx) {
	    int y_off = y_offset + idx;
	    y[y_off] = alpha * x[x_offset + idx] + y[y_off];
	}
    }

    public static void OpXY( int op_len,
			     double[] x, int x_offset,
			     double[] y, int y_offset,
			     IBinaryOp op )
    {
	for (int idx = 0; idx < op_len; ++idx ) {
	    int y_off = y_offset + idx;
	    y[y_off] = op.op(x[x_offset + idx], y[y_off]);
	}
    }

    public static void OpY( int op_len, double[] y, int y_offset, IUnaryOp op )
    {
	for (int idx = 0; idx < op_len; ++idx ) {
	    int y_off = y_offset + idx;
	    y[y_off] = op.op(y[y_off]);
	}
    }
}
