package RScriptReader;
/**
 * function serves to call RScripts from Java currently not used
 */
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class RScriptReader {
	public void read(){
	BufferedReader reader = null;
    Process shell = null;
    try {
        shell = Runtime.getRuntime().exec(new String[] { "/usr/bin/Rscript", "/media/subin/works/subzworks/RLanguage/config/predict.R" });

        reader = new BufferedReader(new InputStreamReader(shell.getInputStream()));
        String line;
        while ((line = reader.readLine()) != null) {
            System.out.println(line);

        }

    } catch (IOException e) {
        e.printStackTrace();
    }
   }
}
