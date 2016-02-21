package netlib_ccm;

import java.util.Arrays;

public final class StridedBuffer
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

    public static StridedBuffer create( int length )
    {
	double[] new_data = new double[length];
	return new StridedBuffer( new_data, 0, length );
    }

    public static void checkRowIndex( int row_count, int row_idx ) throws Exception
    {
	if (row_idx >= row_count)
	    throw new Exception("Attempt to access past end of view");

	if (row_idx < 0)
	    throw new Exception("Attempt to access before beginning of view");
    }

    public static void checkRowIndex( StridedBuffer buffer, int row_idx ) throws Exception
    {
	if (row_idx >= buffer.row_count)
	    throw new Exception("Attempt to access past end of view");

	if (row_idx < 0)
	    throw new Exception("Attempt to access before beginning of view");
    }

    public static int getRowLengthFromRow( int row_count, int column_count
					   , int initial_column_count
					   , int last_column_count
					   , int row )
    {
	if (row == 0)
	    return initial_column_count;
	int last_valid_row = row_count - 1;
	if (row == last_valid_row)
	    return last_column_count;
	return column_count;
    }

    public static int getRowLengthFromRow( StridedBuffer buffer, int row )
    {
	if (row == 0)
	    return buffer.initial_column_count;
	int last_valid_row = buffer.row_count - 1;
	if (row == last_valid_row)
	    return buffer.last_column_count;
	return buffer.column_count;
    }

    public static int getRowLengthFromRowTest(StridedBuffer buffer, int row, int iter_count)
    {
	int initial_column_count = buffer.initial_column_count;
	int row_count = buffer.row_count;
	int last_column_count = buffer.last_column_count;
	int column_count = buffer.column_count;
	int total = 0;
	for (int idx = 0; idx < iter_count; ++idx ) {
	    total += getRowLengthFromRow( row_count, column_count
					  , initial_column_count, last_column_count
					  , row);
	}
	return total;
    }

    public static int getDataLength( int row_count, int column_count
				    , int initial_column_count
				    , int last_column_count )
    {
	return Math.max(0, row_count - 2) * column_count
	    + initial_column_count
	    + last_column_count;
    }

    public static int getDataLength(StridedBuffer buffer)
    {
	return getDataLength( buffer.row_count, buffer.column_count,
			      buffer.initial_column_count, buffer.last_column_count );
    }

    public static void checkDataOffset( int data_length, int data_offset ) throws Exception
    {
	if (data_offset < 0)
	    throw new Exception("Attempt to access before beginning of view");
	if (data_offset >= data_length )
	    throw new Exception("Attempt to access past end of view");
    }

    public static void checkDataOffset(StridedBuffer buffer, int data_offset) throws Exception
    {
	int data_length = getDataLength(buffer);
	checkDataOffset( data_length, data_offset );
    }

    public static int getRowFromDataOffset( int column_count, int initial_column_count
					    , int data_length, int data_offset ) throws Exception
    {
	checkDataOffset( data_length, data_offset );
	if ( data_offset < initial_column_count ) {
	    return 0;
	}
	else {
	    data_offset = data_offset - initial_column_count;
	    return 1 + ( data_offset / column_count );
	}
    }

    public static int getRowFromDataOffset( StridedBuffer buffer,
					    int data_offset ) throws Exception
    {
	return getRowFromDataOffset( buffer.column_count, buffer.initial_column_count
				     , getDataLength(buffer), data_offset );
    }

    public static int getDataOffsetFromRow( int row_count, int column_count
					    , int initial_column_count
					    , int row ) throws Exception
    {
	checkRowIndex( row_count, row );
	if ( row == 0 )
	    return 0;
	return initial_column_count + (column_count * (row - 1));
    }

    public static int getRowLengthFromDataOffset( int row_count, int column_count
						  , int initial_column_count
						  , int last_column_count
						  , int data_offset ) throws Exception
    {
	int data_length = getDataLength( row_count, column_count, initial_column_count
					 , last_column_count );
	checkDataOffset( data_length, data_offset );
	if ( data_offset < initial_column_count )
	    return initial_column_count - data_offset;

	int rest_offset = data_offset - initial_column_count;
	int row_idx = 1 + (rest_offset / column_count);
	int rest_leftover = rest_offset % column_count;
	if ( row_idx == (row_count - 1))
	    return last_column_count - rest_leftover;
	return column_count - rest_leftover;
    }

    public static int getRowLengthFromDataOffset( StridedBuffer buffer
						  , int data_offset ) throws Exception
    {
	return getRowLengthFromDataOffset( buffer.row_count, buffer.column_count
					   , buffer.initial_column_count
					   , buffer.last_column_count
					   , data_offset );
    }

    public static int getTotalOffset( int offset, int row_count, int column_count, int row_stride
				      , int initial_column_count, int last_column_count
				      , int data_offset ) throws Exception
    {
	if ( data_offset == 0 )
	    return offset + Math.max( 0, (column_count
					  - initial_column_count) );

	int data_length = getDataLength( row_count, column_count
					 , initial_column_count, last_column_count );
	checkDataOffset( data_length, data_offset );

	int body_len = column_count * Math.max( 0, row_count - 2);
	if ( data_offset < initial_column_count )
	    return data_offset + offset + (column_count
					   - initial_column_count);
	int data_rest = data_offset - initial_column_count;
	return offset
	    + (row_stride * (1 + (data_rest / column_count)))
	    + (data_rest % column_count);
    }

    public static int getTotalOffset( StridedBuffer buffer, int data_offset ) throws Exception
    {
	return getTotalOffset( buffer.offset, buffer.row_count, buffer.column_count
			       , buffer.row_stride, buffer.initial_column_count
			       , buffer.last_column_count, data_offset );
    }

    public static boolean isDense(StridedBuffer buffer)
    {
	return (1 == buffer.row_count) || (buffer.column_count == buffer.row_stride);
    }

    public static boolean isCompleteDense(StridedBuffer buffer)
    {
	return isDense(buffer)
	    && (buffer.offset == 0)
	    && buffer.data.length == getDataLength(buffer);
    }


    public static int getLongestContiguousLength( int row_count, int column_count
						  , int initial_column_count
						  , int last_column_count
						  , boolean is_dense
						  , int data_length
						  , int data_offset ) throws Exception
    {
	if (is_dense)
	    return Math.max(0, data_length - data_offset);

	return getRowLengthFromDataOffset( row_count, column_count, initial_column_count
					   , last_column_count, data_offset );
    }

    public static int getLongestContiguousLength( StridedBuffer buffer
						  , int data_offset ) throws Exception
    {
	return getLongestContiguousLength( buffer.row_count, buffer.column_count
					   , buffer.initial_column_count
					   , buffer.last_column_count
					   , isDense( buffer )
					   , getDataLength( buffer )
					   , data_offset );
    }

    public StridedBuffer createSubStridedBuffer( int sub_offset, int length ) throws Exception
    {
	if ( length == 0 )
	    return new StridedBuffer( data, offset, 0 );

	int data_len = getDataLength(this);

	if ( length == data_len
	     && sub_offset == 0 )
	    return this;
	checkDataOffset( this, sub_offset );
	checkDataOffset( this, sub_offset + (length - 1));

	int initial_row_len = getLongestContiguousLength( this, sub_offset );
	if ( initial_row_len >= length )
	    return new StridedBuffer( data, getTotalOffset( this, sub_offset ), length );

	int rest_len = length - initial_row_len;
	int num_body_rows  = rest_len / column_count;
	int last_row_len = rest_len % column_count;
	int start_data_offset = column_count - initial_column_count;
	int start_offset = getTotalOffset( this, sub_offset ) - start_data_offset;
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

    public static void unaryOperation( StridedBuffer buffer, IStridedUnaryOp op )
	throws Exception
    {
	int total_op_len = getDataLength(buffer);
	int total_op_offset = 0;
	int offset = buffer.offset;
	int row_count = buffer.row_count;
	int column_count = buffer.column_count;
	int row_stride = buffer.row_stride;
	int initial_column_count = buffer.initial_column_count;
	int last_column_count = buffer.last_column_count;
	boolean is_dense = isDense(buffer);
	while( total_op_offset < total_op_len ) {
	    int current_op_len = getLongestContiguousLength(row_count, column_count
							    , initial_column_count
							    , last_column_count
							    , is_dense
							    , total_op_len
							    , total_op_offset);
	    op.op( buffer.data, getTotalOffset( offset, row_count, column_count, row_stride
						, initial_column_count, last_column_count
						, total_op_offset ), current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public void unaryOperation( IStridedUnaryOp op ) throws Exception
    {
	int total_op_len = getDataLength(this);
	int total_op_offset = 0;
	StridedBuffer buffer = this;
	double[] data = buffer.data;
	while( total_op_offset < total_op_len ) {
	    int current_op_len = getLongestContiguousLength(buffer, total_op_offset);
	    op.op( data, getTotalOffset( buffer, total_op_offset ), current_op_len );
	    total_op_offset += current_op_len;
	}
    }


    public void binaryOperation( StridedBuffer rhs
				 , IStridedBinaryOp op ) throws Exception
    {
	if ( getDataLength(this) == 0
	     || getDataLength(rhs) == 0 )
	    return;
	int total_op_len = getDataLength(this);
	int rhs_len = getDataLength(rhs);
	int total_op_offset = 0;
	while( total_op_offset < total_op_len ) {
	    int rhs_op_offset = total_op_offset % rhs_len;
	    int current_op_len
		= Math.min( getLongestContiguousLength(this, total_op_offset)
			    , rhs.getLongestContiguousLength(rhs, rhs_op_offset) );
	    op.op( data, getTotalOffset(this, total_op_offset)
		   , rhs.data, rhs.getTotalOffset( rhs, rhs_op_offset )
		   , current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public static void binaryOperation( StridedBuffer lhs
					, StridedBuffer rhs
					, IStridedBinaryOp op ) throws Exception
    {
	if ( getDataLength(lhs) == 0
	     || getDataLength(rhs) == 0 )
	    return;
	int total_op_len = getDataLength(lhs);
	int rhs_len = getDataLength(rhs);
	int total_op_offset = 0;

	double[] lhs_data = lhs.data;
	int lhs_offset = lhs.offset;
	int lhs_row_count = lhs.row_count;
	int lhs_column_count = lhs.column_count;
	int lhs_row_stride = lhs.row_stride;
	int lhs_initial_column_count = lhs.initial_column_count;
	int lhs_last_column_count = lhs.last_column_count;
	boolean lhs_is_dense = isDense(lhs);

	double[] rhs_data = rhs.data;
	int rhs_offset = rhs.offset;
	int rhs_row_count = rhs.row_count;
	int rhs_column_count = rhs.column_count;
	int rhs_row_stride = rhs.row_stride;
	int rhs_initial_column_count = rhs.initial_column_count;
	int rhs_last_column_count = rhs.last_column_count;
	boolean rhs_is_dense = isDense(rhs);
	while( total_op_offset < total_op_len ) {
	    int lhs_current_op_len = getLongestContiguousLength(lhs_row_count, lhs_column_count
								, lhs_initial_column_count
								, lhs_last_column_count
								, lhs_is_dense
								, total_op_len
								, total_op_offset);
	    int lhs_total_offset = getTotalOffset( lhs_offset
						   , lhs_row_count, lhs_column_count
						   , lhs_row_stride
						   , lhs_initial_column_count
						   , lhs_last_column_count
						   , total_op_offset );

	    int rhs_op_offset = total_op_offset % rhs_len;
	    int rhs_current_op_len = getLongestContiguousLength(rhs_row_count, rhs_column_count
								, rhs_initial_column_count
								, rhs_last_column_count
								, rhs_is_dense
								, rhs_len
								, rhs_op_offset );
	    int rhs_total_offset = getTotalOffset( rhs_offset
						   , rhs_row_count, rhs_column_count
						   , rhs_row_stride
						   , rhs_initial_column_count
						   , rhs_last_column_count
						   , rhs_op_offset );
	    int current_op_len = Math.min( lhs_current_op_len, rhs_current_op_len );

	    op.op( lhs_data, lhs_total_offset
		   , rhs_data, rhs_total_offset
		   , current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public static void binaryOperationNonLocal( StridedBuffer lhs, StridedBuffer rhs
						, IStridedBinaryOp op ) throws Exception
    {
	if ( getDataLength(lhs) == 0
	     || getDataLength(rhs) == 0 )
	    return;
	int total_op_len = getDataLength(lhs);
	int rhs_len = getDataLength(rhs);
	int total_op_offset = 0;
	double[] lhs_data = lhs.data;
	double[] rhs_data = rhs.data;
	while( total_op_offset < total_op_len ) {
	    int rhs_op_offset = total_op_offset % rhs_len;
	    int current_op_len
		= Math.min( getLongestContiguousLength(lhs, total_op_offset)
			    , rhs.getLongestContiguousLength(rhs, rhs_op_offset) );
	    op.op( lhs_data, getTotalOffset(lhs, total_op_offset)
		   , rhs_data, rhs.getTotalOffset( rhs, rhs_op_offset )
		   , current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public static void ternaryOperation( StridedBuffer lhs, StridedBuffer middle
					 , StridedBuffer rhs
					 , IStridedTernaryOp op ) throws Exception
    {
	int total_op_len = getDataLength(lhs);
	int middle_len = getDataLength(middle);
	int rhs_len = getDataLength(rhs);
	if ( total_op_len == 0
	     || middle_len == 0
	     || rhs_len == 0 )
	    return;

	int total_op_offset = 0;

	double[] lhs_data = lhs.data;
	int lhs_offset = lhs.offset;
	int lhs_row_count = lhs.row_count;
	int lhs_column_count = lhs.column_count;
	int lhs_row_stride = lhs.row_stride;
	int lhs_initial_column_count = lhs.initial_column_count;
	int lhs_last_column_count = lhs.last_column_count;
	boolean lhs_is_dense = isDense(lhs);

	double[] rhs_data = rhs.data;
	int rhs_offset = rhs.offset;
	int rhs_row_count = rhs.row_count;
	int rhs_column_count = rhs.column_count;
	int rhs_row_stride = rhs.row_stride;
	int rhs_initial_column_count = rhs.initial_column_count;
	int rhs_last_column_count = rhs.last_column_count;
	boolean rhs_is_dense = isDense(rhs);


	double[] middle_data = middle.data;
	int middle_offset = middle.offset;
	int middle_row_count = middle.row_count;
	int middle_column_count = middle.column_count;
	int middle_row_stride = middle.row_stride;
	int middle_initial_column_count = middle.initial_column_count;
	int middle_last_column_count = middle.last_column_count;
	boolean middle_is_dense = isDense(middle);

	while( total_op_offset < total_op_len ) {
	    int lhs_current_op_len = getLongestContiguousLength(lhs_row_count, lhs_column_count
								, lhs_initial_column_count
								, lhs_last_column_count
								, lhs_is_dense
								, total_op_len
								, total_op_offset);
	    int lhs_total_offset = getTotalOffset( lhs_offset
						   , lhs_row_count, lhs_column_count
						   , lhs_row_stride
						   , lhs_initial_column_count
						   , lhs_last_column_count
						   , total_op_offset );

	    int rhs_op_offset = total_op_offset % rhs_len;
	    int rhs_current_op_len = getLongestContiguousLength(rhs_row_count, rhs_column_count
								, rhs_initial_column_count
								, rhs_last_column_count
								, rhs_is_dense
								, rhs_len
								, rhs_op_offset );

	    int rhs_total_offset = getTotalOffset( rhs_offset
						   , rhs_row_count, rhs_column_count
						   , rhs_row_stride
						   , rhs_initial_column_count
						   , rhs_last_column_count
						   , rhs_op_offset );

	    int middle_op_offset = total_op_offset % middle_len;
	    int middle_current_op_len = getLongestContiguousLength(middle_row_count
								   , middle_column_count
								   , middle_initial_column_count
								   , middle_last_column_count
								   , middle_is_dense
								   , middle_len
								   , middle_op_offset );

	    int middle_total_offset = getTotalOffset( middle_offset
						      , middle_row_count, middle_column_count
						      , middle_row_stride
						      , middle_initial_column_count
						      , middle_last_column_count
						      , middle_op_offset );


	    int current_op_len = Math.min( lhs_current_op_len
					   , Math.min( rhs_current_op_len
						       , middle_current_op_len ) );
	    op.op( lhs_data, lhs_total_offset
		   , middle_data, middle_total_offset
		   , rhs_data, rhs_total_offset
		   , current_op_len );
	    total_op_offset += current_op_len;
	}
    }

    public static void unaryOperationOp( StridedBuffer lhs
					 , IUnaryOp unary_op ) throws Exception
    {
	unaryOperation( lhs, new IStridedUnaryOp() {
		public void op( double[] data, int offset, int len ) {
		    Ops.OpY( len, data, offset, unary_op );
		}
	    } );
    }

    public static void binaryOperationOp( StridedBuffer lhs, StridedBuffer rhs,
					  IBinaryOp bin_op ) throws Exception
    {
	int lhs_len = getDataLength(lhs);
	int rhs_len = getDataLength(rhs);
	if ( rhs_len == 0 || lhs_len == 0 )
	    return;
	if ( rhs_len == 1 ) {
	    double rhs_val = rhs.data[getTotalOffset(rhs, 0)];
	    IUnaryOp condensed_op = new IUnaryOp() {
		    public double op( double lhs_val ) {
			return bin_op.op( lhs_val, rhs_val );
		    }
		};
	    unaryOperation( lhs, new IStridedUnaryOp() {
		    public void op( double[] data, int offset, int len ) {
			Ops.OpY( len, data, offset, condensed_op );
		    }
		} );
	} else {
	    binaryOperation( lhs, rhs, new IStridedBinaryOp() {
		    public void op( double[] lhs_data, int lhs_offset
				    , double[] rhs_data, int rhs_offset
				    , int op_len ) {
			Ops.OpXY( op_len, rhs_data, rhs_offset,
				  lhs_data, lhs_offset, bin_op );
		    }
		} );
	}
    }

    public void set( double value ) throws Exception
    {
	try {
	    unaryOperation( this, new IStridedUnaryOp() {
		    public void op( double[] lhs, int lhs_offset, int op_len ) {
			Arrays.fill( lhs, lhs_offset, lhs_offset + op_len, value );
		    }
		} );
	}
	//Not going to happen
	catch( Exception e ) {
	}
    }

    //Assign other to this
    public void assign( StridedBuffer other ) throws Exception
    {
	if ( getDataLength( other ) == 1 ) {
	    double val = other.data[getTotalOffset(other, 0)];
	    set( val );
	}
	else {
	    binaryOperation( other, new IStridedBinaryOp() {
		    public void op( double[] lhs, int lhs_offset
				    , double[] rhs, int rhs_offset, int op_len ) {
			System.arraycopy( rhs, rhs_offset, lhs, lhs_offset, op_len );
		    }
		} );
	}
    }

    //Create new *dense* strided buffer with same data as this one.
    public StridedBuffer clone()
    {
	StridedBuffer other = create(getDataLength(this));
	try
	{
	    other.assign( this );
	    return other;
	}
	//will never happen.
	catch( Exception e ) {}
	//But if it does I do want to know about it.
	return null;
    }

    public static boolean areOverlapping( StridedBuffer lhs, StridedBuffer rhs )
	throws Exception
    {
	if ( lhs.data == rhs.data
	     && getDataLength( lhs ) > 0
	     && getDataLength( rhs ) > 0 ) {
	    int lhs_start = getTotalOffset( lhs, 0 );
	    int rhs_start = getTotalOffset( rhs, 0 );
	    int lhs_end = getTotalOffset( lhs, ( getDataLength( lhs ) - 1 ) );
	    int rhs_end = getTotalOffset( rhs, ( getDataLength( rhs ) - 1 ) );
	    if ( ( lhs_start >= rhs_start && lhs_start <= rhs_end )
		 || ( rhs_start >= lhs_start && rhs_start <= lhs_end ) )
		return true;
	}
	return false;
    }

    public static void MulTest( StridedBuffer lhs, StridedBuffer rhs, int iter_count )
    {
	double[] rhs_data = rhs.data;
	double[] lhs_data = lhs.data;
	int op_len = rhs_data.length;
	for ( int idx = 0; idx < iter_count; ++idx ) {
	    Ops.AY( op_len, rhs_data, 0, lhs_data, 0 );
	}
    }

}
