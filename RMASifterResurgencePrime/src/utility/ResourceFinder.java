package utility;

import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.FileVisitor;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
/**
 * Tjis class is used to locate the NCBI Map and NCBI tree but it is way faster 
 * to provide the file Location to MaltExtract
 * @author huebler
 *
 */
public class ResourceFinder {
	public String getPath(String pattern){
		ArrayList<Path> files = new ArrayList<Path>();
		Path startDir = Paths.get("/");
		FileSystem fs = FileSystems.getDefault();
		String currentFile = "";
		final PathMatcher matcher = fs.getPathMatcher("glob:" + pattern);

		FileVisitor<Path> matcherVisitor = new SimpleFileVisitor<Path>() {
			
		    // Invoke the pattern matching
	        // method on each directory.
	        @Override
	        public FileVisitResult preVisitDirectory(Path dir,
	                BasicFileAttributes attrs) {
	        	// System.out.print(String.format("Found matched file: '%s'.%n", dir));
	            return FileVisitResult.CONTINUE;
	        }

	        @Override
	        public FileVisitResult visitFileFailed(Path file,
	                IOException exc) {
	           // System.err.println(exc);
	            return FileVisitResult.CONTINUE;
	        }
		    @Override
		    public FileVisitResult visitFile(Path file, BasicFileAttributes attribs) {
		        Path name = file.getFileName();
		        if (matcher.matches(name)) {
		           files.add(file);
		        }
		        return FileVisitResult.CONTINUE;
		    }
		};
		try {
			Files.walkFileTree(startDir, matcherVisitor);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		long mostRecent = 0;
		for(Path path : files){
			long fileTime = path.toFile().lastModified();
			if(fileTime>=mostRecent){
			mostRecent =	fileTime;
			currentFile = path.toString();
			}
		}
		return currentFile;
	}
	
}
