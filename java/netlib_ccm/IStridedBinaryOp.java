package netlib_ccm;


public interface IStridedBinaryOp
{
    public void op(double[] lhs_data, int lhs_offset, double[] rhs_Data, int rhs_offset, int op_len);
}
