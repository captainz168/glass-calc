# Glass Calc — โปรแกรมคำนวณกระจก

เครื่องมือภายในบริษัท Kama Joint Venture สำหรับคำนวณความหนากระจกรับแรงลมและความกว้างครีบกระจก (Rib)
เป็นไฟล์ HTML ไฟล์เดียว เปิดใน browser ได้เลย ใช้งาน offline ได้ ไม่ต้องติดตั้ง

| ไฟล์ | คืออะไร |
|---|---|
| `GL300_new.html` | คำนวณตาม **ASTM E1300** (ลม 3 วินาที, Pb = 8/1000): LR = NFL × GTF × LS |
| `gl300-legacy.html` | สูตรและตารางเดียวกับโปรแกรม **GL300 v2.0** เดิม (แกะจากตัวโปรแกรม) แก้ bug ช่วง b/a 1.8–2.0 แล้ว |
| `legacy/GL300/` | โปรแกรม GL300 v2.0 ตัวเดิม (VB4, 1997) เก็บไว้อ้างอิง |
| `src/index.src.html` | ต้นฉบับของ `GL300_new.html` (ยังไม่ใส่ตารางข้อมูล) |
| `tools/` | สคริปต์ Python ที่ใช้สร้างตาราง NFL และ `tools/data/` ตารางที่สร้างแล้ว |

ทั้งสองไฟล์รองรับภาษาไทย / English, การยึด 4 / 3 / 2 ด้าน, กระจก Annealed / Heat Strengthened / Tempered, Laminated, Insulated, ครีบกระจก, พิมพ์ใบคำนวณ และบันทึกงานเป็น .json

## GL300_new.html — ที่มาของค่า NFL

| การยึด | วิธี |
|---|---|
| 4 ด้าน กระจกชั้นเดียว | วิธีวิเคราะห์ของ E1300 (§6.3): แผ่น nonlinear (von Kármán) + failure prediction model (m = 7, k = 2.86×10⁻⁵³ N⁻⁷m¹²) ตรวจทับกราฟ Fig. A1.1–A1.14 แล้ว |
| 4 ด้าน Laminated, 3 ด้าน | ค่าจากกราฟ E1300 Fig. A1.15–A1.26, A1.29–A1.42 (อ่านเส้นด้วยโปรแกรมแล้วตรวจทับกราฟต้นฉบับ) |
| 2 ด้าน | Fig. A1.27: NFL = 21,600·t²/L² (kPa, mm); Laminated จาก Fig. A1.43 |

อ้างอิงหลักเป็น ASTM E1300-16 (ตาราง GTF ของ Insulated จาก E1300-04) **ยังไม่ได้เทียบกับ E1300-24**
ผลเป็นเครื่องมือช่วยออกแบบ งานจริงควรให้วิศวกรตรวจรับรอง

## สร้าง GL300_new.html ใหม่ (หลังแก้ `src/index.src.html`)

```bash
pip install numpy
python3 tools/build.py tools/data GL300_new.html
```

## เอกสารอ้างอิง
1. ASTM E1300-16, *Standard Practice for Determining Load Resistance of Glass in Buildings* — §6.2, §6.3, Tables 1 & 4, Annex A1, Annex A2, Appendix X9
2. ASTM E1300-04, Table 2 (GTF insulating glass) และ Appendix X3 (load share)
3. Beason, W.L. & Morgan, J.R. (1984), "Glass Failure Prediction Model", *J. Struct. Eng.*, ASCE, 110(2), 197–212
