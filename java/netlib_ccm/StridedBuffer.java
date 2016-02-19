package netlib_ccm;


public class StridedBuffer
{
    public final double[] data;
    public final int offset;
    public final int row_count;
    public final int column_count;
    public final int row_stride;
    public final int initial_column_count;
    public final int last_column_count;
    public StridedBuffer( double[] data, int offset
			  , int row_count, int column_count
			  , int row_stride, int initial_column_count
			  , int last_column_count )
    {
	this.data = data;
	this.offset = offset;
	this.row_count = row_count;
	this.column_count = column_count;
	this.row_stride = row_stride;
	this.initial_column_count = initial_column_count;
	this.last_column_count = last_column_count;
    }

    public StridedBuffer( double[] data, int offset, int length )
    {
	this.data = data;
	this.offset = offset;
	this.row_count = 1;
	this.column_count = length;
	this.row_stride = length;
	this.initial_column_count = length;
	this.last_column_count = 0;
    }

    public void checkRowIndex( int row_idx ) throws Exception
    {
	if (row_idx >= row_count)
	    throw new Exception("Attempt to access past end of view");

	if (row_idx < 0)
	    throw new Exception("Attempt to access before beginning of view");
    }

    public int getRowLengthFromRow( int row )
    {
	int last_valid_row = row_count - 1;
	if (row == 0)
	    return initial_column_count;
	if (row == last_valid_row)
	    return last_column_count;
	return column_count;
    }

    public int getDataLength()
    {
	return Math.max(0, row_count - 2) * column_count
	    + initial_column_count
	    + last_column_count;
    }

    public void checkDataOffset(int data_offset) throws Exception
    {
	if (data_offset < 0)
	    throw new Exception("Attempt to access before beginning of view");
	if (data_offset >= getDataLength())
	    throw new Exception("Attempt to access past end of view");
    }

    public int getRowFromDataOffset( int data_offset ) throws Exception
    {
	checkDataOffset( data_offset );
	if ( data_offset < initial_column_count ) {
	    return 0;
	}
	else {
	    data_offset = data_offset - initial_column_count;
	    return 1 + ( data_offset / column_count );
	}
    }

    public int getDataOffsetFromRow( int row ) throws Exception
    {
	checkRowIndex( row );
	if ( row == 0 )
	    return 0;
	return initial_column_count + (column_count * (row - 1));
    }

    public int getRowLengthFromDataOffset( int data_offset ) throws Exception
    {
	checkDataOffset( data_offset );
	if ( data_offset < initial_column_count )
	    return initial_column_count - data_offset;

	int rest_offset = data_offset - initial_column_count;
	int row_idx = 1 + (rest_offset / column_count);
	int rest_leftover = rest_offset % column_count;
	if ( row_idx == (row_count - 1))
	    return last_column_count - rest_leftover;
	return column_count - rest_leftover;
    }

    public int getTotalOffset( int data_offset ) throws Exception
    {
	if ( data_offset == 0 )
	    return offset + Math.max( 0, (column_count - initial_column_count) );
	checkDataOffset( data_offset );
	int data_len = getDataLength();
	int body_len = column_count * Math.max( 0, row_count - 2);
	if ( data_offset < initial_column_count )
	    return data_offset + offset + (column_count - initial_column_count);
	int data_rest = data_offset - initial_column_count;
	return offset
	    + (row_stride * (1 + (data_rest / column_count)))
	    + (data_rest % column_count);
    }

    public bool isDense()
    {
	return (1 == row_count) || (column_count == row_stride);
    }

    public bool isCompleteDense()
    {
	return isDense()
	    && (offset == 0)
	    && (.length data) == getViewDataLength();
    }

    public StridedBuffer createSubStridedBuffer( int sub_offset, int length )
    {
	if ( length == 0 )
	    return new StridedBuffer( data, 0, 0 );

	if ( length == getDataLength()
	     && sub_offset == 0 )
	    return this;
	checkDataOffset( sub_offset );
	checkDataOffset( sub_offset + (length - 1));
	int data_len = getDataLength();
	int first_row_len = 0;
	if ( sub_offset < initial_column_count ) {
	    first_row_len = initial_column_count - sub_offset;
	}
	else {
	    first_row_len = column_count - ((sub_offset - initial_column_count) % column_count);
	}
	int first_row_offset = column_count - first_row_len;

	if ( length < first_row_len )
	    return new StridedBuffer( data, getTotalOffset( sub_offset ), length );

	int row_idx = getRowFromDataOffset(sub_offset);
	int last_row_len = ((first_row_offset + length) % column_count);
	int ary_offset (offset + (row_stride * row-idx));
	int body_num_rows = ((length - (first_row_len + last_row_len)) / column_count);

	int ini_row = 0;
	if ( first_row_len != 0 ) ini_row = 1;

	int end_row = 1;
	return new StridedBuffer( data, ary_offset
				  , (body_num_rows + ini_row + end_row)
				  , column_count
				  , row_stride
				  , first_row_len
				  , last_row_len );
    }

    public void stridedOperation( StridedBuffer rhs, IStridedBinaryOp op ) throws Exception
    {
	if ( getDataLength() == 0
	     || rhs.getDataLength() == 0 )
	    return;
	int num_ops = getDataLength / rhs.getDataLength();
	int op_len = rhs.getDataLength();
	bool are_both_dense = isDense() && rhs.isDense();
	bool rhs_dense = rhs.isDense();
	for ( int op_idx = 0; op_idx < num_ops; ++op_idx ) {
	    int lhs_offset = op_len * op_idx;

	    if ( are_both_dense ) {
		op.op( data, getTotalOffset(lhs_offset), rhs.data, rhs.getTotalOffset( 0 ), op_len );
	    }
	    else {
		if ( rhs_dense ) {
		    int rhs_offset = 0;
		    //Iterate lhs rows because rhs is dense
		    for ( int lhs_row = getRowFromDataOffset(lhs_offset);
			  lhs_row < row_count && rhs_offset < op_len; ++lhs_row ) {
			int lhs_row_len = getRowLengthFromDataOffset(lhs_offset + rhs_offset);
			int max_possible = op_len - rhs_offset;
			int copy_len = Math.min(lhs_row_len, max_possible);
			op( data, getTotalOffset( lhs_offset + rhs_offset ),
			    rhs.data, rhs.getTotalOffset( rhs_offset ),
			    copy_len );
			rhs_offset += copy_len;
		    }
		}
		else {
		    //Harder case, iterate rhs rows because both may be strided
		    int rhs_offset = 0;
		    while( rhs_offset < op_len ) {
			int rhs_row_len = getRowLengthFromDataOffset( rhs_offset );
			int lhs_row_len = getRowLengthFromDataOffset( lhs_offset + rhs_offset );
			int copy_len = Math.min( rhs_row_len, lsh_row_len );
			op( data, getTotalOffset( lhs_offset + rhs_offset )
			    , rhs.data, rhs.getTotalOffset( rhs_offset )
			    , copy_len );
			rhs_offset += copy_len;
		    }
		}
	    }
	}
    }
}
