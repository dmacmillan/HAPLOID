$(document).ready(function() {

// Loading samples 
// **************************************************************
$('#load_sample').click(function() {
	sequences = 'PAT1	A*01:01	A*02:01	B*07:02	B*08:01	C*01:02	C*01:02	CAACAAAAAAAA\n' +
                'PAT2	A*02:01	A*03:01	B*13:02	B*57:01	C*01:02	C*06:02	AAAAAAGAAAAA\n' +
                'PAT3	A*11:01	A*23:01	B*14:01	B*57:01	C*01:02	C*06:02	AAAAAAGAAAAA\n' +
                'PAT4	A*01:01	A*02:01	B*07:02	B*15:01	C*01:02	C*01:02	AAACAAAAAAAA\n' +
                'PAT5	A*03:01	A*11:01	B*08:01	B*57:01	C*01:02	C*06:02	AAAAAAGAAAAA\n' +
                'PAT6	A*01:01	A*23:01	B*13:02	B*14:01	C*01:02	C*01:02	AAACAAAAAAAA\n' +
                'PAT7	A*02:01	A*03:01	B*15:01	B*57:01	C*01:02	C*06:02	AAAAAAGAAAAA\n' +
                'PAT8	A*01:01	A*11:01	B*07:02	B*08:01	C*01:02	C*01:02	AAACAAAAAAAA\n' +
                'PAT9	A*02:01	A*23:01	B*13:02	B*57:01	C*01:02	C*06:02	AAAAAAGAAAAA\n' +
                'PAT10	A*03:01	A*11:01	B*14:01	B*15:01	C*01:02	C*06:02	AAAAAAAKAAAA'
	$('#sequences').val(sequences);
    $('#hlacountfilter').val(1);
    $('#aacountfilter').val(1);
    $('#maxp').val(0.05);
    $('#maxq').val(0.2);
});
// **************************************************************

// Clearing
// **************************************************************
$('#clear').click(function() {
	$('#sequences').val('');
	$('#hlacountfilter').val('');
	$('#aacountfilter').val('');
	$('#maxp').val('');
	$('#maxq').val('');
});
// **************************************************************

});
