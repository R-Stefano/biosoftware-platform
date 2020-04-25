var restify = require('restify');
const path=require('path')

const server = restify.createServer();
const port=process.env.port || 8080
//Serve frontend as static
server.get('/*', restify.plugins.serveStatic({
	directory: path.join(__dirname, 'webapp/dist/frontend/'),
	default: 'index.html'
}));

server.listen(port, () => {
	console.log('Service server started. Listening on port '+ port);
});