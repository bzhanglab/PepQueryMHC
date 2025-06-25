package progistar.scan.data;

import java.util.ArrayList;
import java.util.Hashtable;

public class LocTable {
	public Hashtable<String, Hashtable<String, LocationInformation>> table = new Hashtable<String, Hashtable<String, LocationInformation>>(); 

	public boolean putLocation (LocationInformation lInfo) {
		Hashtable<String, LocationInformation> gTable = table.get(lInfo.inputSequence);
		if(gTable == null) {
			gTable = new Hashtable<String, LocationInformation>();
			table.put(lInfo.inputSequence, gTable);
		}
		
		LocationInformation slInfo = gTable.get(lInfo.getKey());
		if(slInfo == null) {
			gTable.put(lInfo.getKey(), lInfo);
			return true;
		} else {
			lInfo.readCounts.forEach((barcodeId, value)->{
				Long val = slInfo.readCounts.get(barcodeId);
				if(val == null) {
					val = 0L;
				}
				slInfo.readCounts.put(barcodeId, (value + val));
			});
		}
		
		return false;
	}
	
	public ArrayList<LocationInformation> getLocations (String inputSequence) {
		ArrayList<LocationInformation> locations = new ArrayList<LocationInformation>();
		Hashtable<String, LocationInformation> gTable = table.get(inputSequence);
		if(gTable != null) {
			gTable.forEach((key, info)->{
				locations.add(info);
			});
		}
		
		// if there is no matched reads
		if(locations.size() == 0) {
			LocationInformation nullLocation = new LocationInformation();
			nullLocation.readCounts.put(Constants.DEFAULT_BARCODE_ID, 0L);
			locations.add(nullLocation);
		}
		
		return locations;
	}
}
