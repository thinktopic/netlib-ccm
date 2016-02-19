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

    public boolean isDense()
    {
	return (1 == row_count) || (column_count == row_stride);
    }

    public boolean isCompleteDense()
    {
	return isDense()
	    && (offset == 0)
	    && data.length == getDataLength();
    }


    public int getLongestContiguousLength( int data_offset ) throws Exception
    {
	if (isDense())
	    return Math.max(0, getDataLength() - data_offset);

	return getRowLengthFromDataOffset( data_offset );
    }

    public StridedBuffer createSubStridedBuffer( int sub_offset, int length ) throws Exception
    {
	if ( length == 0 )
	    return new StridedBuffer( data, offset, 0 );

	int data_len = getDataLength();

	if ( length == data_len
	     && sub_offset == 0 )
	    return this;
	checkDataOffset( sub_offset );
	checkDataOffset( sub_offset + (length - 1));

	int initial_row_len = getLongestContiguousLength( sub_offset );
	if ( initial_row_len >= length )
	    return new StridedBuffer( data, getTotalOffset( sub_offset ), length );

	int rest_len = length - initial_row_len;
	int num_body_rows  = rest_len / column_count;
	int last_row_len = rest_len % column_count;
	int start_data_offset = column_count - initial_column_count;
	int start_offset = getTotalOffset( sub_offset ) - start_data_offset;
	int num_rows = num_body_rows;
	if (initial_row_len > 0)
	    num_rows++;
	if (last_row_len > 0)
	    num_rows++;

	return new StridedBuffer( data, start_offset, num_rows
				  , column_count
				  , row_stride
				  , initial_row_len
				  , last_row_len );
    }

    public void unaryOperation( IStridedUnaryOp op ) throws Exception
    {
	int total_op_len = getDataLength();
	int total_op_offset = 0;
	while( total_op_offset < total_op_len ) {
	    int current_op_len = getLongestContiguousLength(total_op_offset);
	    op.op( data, getTotalOffset( total_op_offset ), current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public void binaryOperation( StridedBuffer rhs, IStridedBinaryOp op ) throws Exception
    {
	if ( getDataLength() == 0
	     || rhs.getDataLength() == 0 )
	    return;
	int total_op_len = getDataLength();
	int rhs_len = rhs.getDataLength();
	int total_op_offset = 0;
	while( total_op_offset < total_op_len ) {
	    int rhs_op_offset = total_op_offset % rhs_len;
	    int current_op_len = Math.min( getLongestContiguousLength(total_op_offset)
					   , rhs.getLongestContiguousLength(rhs_op_offset) );
	    op.op( data, getTotalOffset(total_op_offset)
		   , rhs.data, rhs.getTotalOffset( rhs_op_offset )
		   , current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public void ternaryOperation( StridedBuffer middle, StridedBuffer rhs, IStridedTernaryOp op ) throws Exception
    {
	if ( getDataLength() == 0
	     || middle.getDataLength() == 0
	     || rhs.getDataLength() == 0 )
	    return;
	int total_op_len = getDataLength();
	int middle_len = middle.getDataLength();
	int rhs_len = rhs.getDataLength();
	int total_op_offset = 0;
	while( total_op_offset < total_op_len ) {
	    int middle_op_offset = total_op_offset % middle_len;
	    int rhs_op_offset = total_op_offset % rhs_len;
	    int current_op_len = Math.min( getLongestContiguousLength(total_op_offset)
					   , Math.min( middle.getLongestContiguousLength( middle_op_offset )
						       , rhs.getLongestContiguousLength( rhs_op_offset ) ) );
	    op.op( data, getTotalOffset( total_op_offset )
		   , middle.data, middle.getTotalOffset( middle_op_offset )
		   , rhs.data, rhs.getTotalOffset( rhs_op_offset )
		   , current_op_len );
	    total_op_offset += current_op_len;
	}
    }
}
