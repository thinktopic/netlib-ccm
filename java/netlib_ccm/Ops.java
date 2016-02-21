package netlib_ccm;


public class Ops
{
    public static void alphaAXPlusBetaY( int op_len, double alpha,
					 double[] a, int a_offset,
					 double[] x, int x_offset,
					 double beta,
					 double[] y, int y_offset)
    {
	boolean zero_offsets = a_offset == 0 && x_offset == 0 && y_offset == 0;
	if ( beta == 1.0 ) {
	    if ( zero_offsets ) {
		for (int idx = 0; idx < op_len; ++idx) {
		    y[idx] += alpha * x[idx] * a[idx];
		}
	    }
	    else {
		for (int idx = 0; idx < op_len; ++idx) {
		    int y_off = y_offset + idx;
		    y[y_off] += alpha * x[x_offset + idx] * a[a_offset + idx];
		}
	    }
	}
	else {
	    if ( alpha == 1.0 ) {
		if ( zero_offsets ) {
		    for (int idx = 0; idx < op_len; ++idx) {
			y[idx] = x[idx] * a[idx] + beta * y[idx];
		    }
		}
		else {
		    for (int idx = 0; idx < op_len; ++idx) {
			int y_off = y_offset + idx;
			y[y_off] = x[x_offset + idx] * a[a_offset + idx] + beta * y[y_off];
		    }
		}
	    }
	    else {
		if ( zero_offsets ) {
		    for (int idx = 0; idx < op_len; ++idx) {
			y[idx] = alpha * x[idx] * a[idx] + beta * y[idx];
		    }
		}
		else {
		    for (int idx = 0; idx < op_len; ++idx) {
			int y_off = y_offset + idx;
			y[y_off] = alpha * x[x_offset + idx] * a[a_offset + idx]
			    + beta * y[y_off];
		    }
		}
	    }
	}
    }

    public static void AY( int op_len,
			   double[] a, int a_offset,
			   double[] y, int y_offset )
    {
	if ( a_offset == 0 && y_offset == 0 ) {
	    for ( int idx = 0; idx < op_len; ++idx )
		y[idx] *= a[idx];
	}
	else {
	    for (int idx = 0; idx < op_len; ++idx) {
		y[y_offset + idx] *= a[a_offset + idx];
	    }
	}
    }

    public static void alphaAY( int op_len, double alpha,
				double[] a, int a_offset,
				double[] y, int y_offset )
    {
	if ( alpha == 1.0 ) {
	    AY(op_len, a, a_offset, y, y_offset);
	}
	else {
	    if ( a_offset == 0 && y_offset == 0 ) {
		for ( int idx = 0; idx < op_len; ++idx )
		    y[idx] *= alpha * a[idx];
	    }
	    for (int idx = 0; idx < op_len; ++idx) {
		int y_off = y_offset + idx;
		y[y_off] *= alpha * a[a_offset + idx];
	    }
	}
    }

    public static void axpy( int op_len, double alpha,
			     double[] x, int x_offset,
			     double[] y, int y_offset )
    {
	if ( x_offset == 0 && y_offset == 0 ) {
	    for (int idx = 0; idx < op_len; ++idx) {
		y[idx] += alpha * x[idx];
	    }
	}
	else {
	    for (int idx = 0; idx < op_len; ++idx) {
		int y_off = y_offset + idx;
		y[y_off] += alpha * x[x_offset + idx];
	    }
	}
    }

    public static void DivXY( int op_len,
			      double[] x, int x_offset,
			      double[] y, int y_offset )

    {
	if ( y_offset == 0 && x_offset == 0 ) {
	    for (int idx = 0; idx < op_len; ++idx ) {
		y[idx] /= x[idx];
	    }
	}
	else {
	    for (int idx = 0; idx < op_len; ++idx ) {
		y[y_offset + idx] /= x[x_offset + idx];
	    }
	}
    }

    public static void OpXY( int op_len,
			     double[] x, int x_offset,
			     double[] y, int y_offset,
			     IBinaryOp op )
    {
	if ( x_offset == 0 && y_offset == 0 ) {
	    for (int idx = 0; idx < op_len; ++idx ) {
		y[idx] = op.op(x[idx], y[idx]);
	    }
	}
	else {
	    for (int idx = 0; idx < op_len; ++idx ) {
		int y_off = y_offset + idx;
		y[y_off] = op.op(x[x_offset + idx], y[y_off]);
	    }
	}
    }

    public static void OpY( int op_len, double[] y, int y_offset, IUnaryOp op )
    {
	if ( y_offset == 0 ) {
	    for (int idx = 0; idx < op_len; ++idx ) {
		y[idx] = op.op(y[idx]);
	    }
	}
	else {
	    for (int idx = 0; idx < op_len; ++idx ) {
		int y_off = y_offset + idx;
		y[y_off] = op.op(y[y_off]);
	    }
	}
    }
}
