package progistar.scan.data;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ParseCigar {
	private static final Pattern EACH_CIGAR_REGEX = Pattern.compile("([0-9]+)([MINDSHPX=])");
	
	public static void parseCigar (String cigar) {
		Matcher matcher = EACH_CIGAR_REGEX.matcher(cigar);
		while (matcher.find()) {
			int markerSize = Integer.parseInt(matcher.group(1));
			char operation = matcher.group(2).charAt(0);
	      
	    }
	}
}
