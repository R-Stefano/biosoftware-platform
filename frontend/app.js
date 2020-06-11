var express = require('express');
const app = express()
const path=require('path')

// Point static path to dist
app.use(express.static(path.join(__dirname, 'webapp/dist/frontend')));

//Serve frontend as static
app.get('*', (req,res) =>{
    res.sendFile(path.join(__dirname,'webapp/dist/frontend/index.html'));
});

app.listen(process.env.port || 8080, () => {
	console.log('Service server started. Listening on port '+(process.env.port || 8080));
});