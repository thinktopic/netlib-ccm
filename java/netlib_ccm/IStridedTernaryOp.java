package netlib_ccm;


public interface IStridedTernaryOp
{
    public void op( double[] lhs_data, int lhs_offset
		    , double[] middle_data, int middle_offset
		    , double[] rhs_data, int rhs_offset
		    , int op_len );
}
