QString PColorsTab::load_doc()
{
	return
	"<html><body>"
	"<br>"
	"<center>"
	"<h2> Paricle Color Codes</h2>"
	"</center>"
	"<br><br>"
	"<table bordersize='0'>"
	"<tr><td>Particle</td>            <td>Color</td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td>photons:</td>            <td><font color='blue'>   blue   </font></td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td>e<sup>-</sup>:  </td>               <td><font color='cyan'>   cyan   </font></td></tr>"
	"<tr><td>e<sup>+</sup>:  </td>               <td><font color='#00CCCC'>   dark cyan   </font></td></tr>"
	"<tr><td>&nu;<sub>e</sub>:</td>   <td><font color='#99CCFF'> azure </font></td></tr>"
	"<tr><td>anti&nu;<sub>e</sub>:</td>   <td><font color='#0080FF'> dark azure </font></td></tr>"

	"<tr><td>&mu;<sup>+</sup>:</td>   <td><font color='#CC99FF'> violet </font></td></tr>"
	"<tr><td>&mu;<sup>-</sup>:</td>   <td><font color='#7F00FF'> dark violet </font></td></tr>"
	"<tr><td>&nu;<sub>&mu;</sub>:</td>   <td><font color='#FFCCFF'> pink </font></td></tr>"
	"<tr><td>anti&nu;<sub>&mu;</sub>:</td>   <td><font color='#CC00CC'> dark pink </font></td></tr>"
	"<tr><td><br></td><td></td></tr>"

	"<tr><td>protons: </td>           <td><font color='orange'> orange </font></td></tr>"
	"<tr><td>neutrons:</td>           <td><font color='black'>  black </font></td></tr>"
	"<tr><td><br></td><td></td></tr>"

	"<tr><td>&pi;<sup>+</sup>:</td>   <td><font color='magenta'> magenta </font></td></tr>"
	"<tr><td>&pi;<sup>-</sup>:</td>   <td><font color='yellow'>  yellow  </font></td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td colspan='2'>all other particle colors are based on their charges:</td></tr>"
	"<tr><td><br></td><td></td></tr>"
	"<tr><td>q=+1:</td>               <td><font color='red'>   red   </font></td></tr>"
	"<tr><td>q= 0:</td>               <td><font color='black'> white </font></td></tr>"
	"<tr><td>q=-1:</td>               <td><font color='green'> green </font></td></tr>"
	"</table>"
	" </body></html>";
}
